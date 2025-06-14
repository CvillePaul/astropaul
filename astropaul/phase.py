from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass, replace

from astropy.time import Time
import astropy.units as u
import pandas as pd


@dataclass
class Ephemeris:
    """Holds all parameters necessary to calculate eclipse events"""

    system: str
    component: str
    t0: float
    period: float
    duration: float
    # eccentricity: float

    def from_dataframe_row(row: pd.DataFrame) -> "Ephemeris":
        return Ephemeris(*row[["System", "Member", "T0", "Period", "Duration"]])


@dataclass
class PhaseEvent:
    """Meant to encode the details of one particular occurrence in the lifetime of a bound system.
    Examples include ingress, egress, mid eclipse, but could be anything defined by a PhaseEventDef"""

    time: Time  # a value of None means the event occurred sometime before the current time segment
    system: str
    member: str
    orbit: int
    type: str
    phase: float

    def copy(self) -> "PhaseEvent":
        return PhaseEvent(self.time, self.system, self.member, self.orbit, self.type, self.phase)

    def to_list(self) -> list:
        return [self.system, self.member, self.orbit, self.time, self.phase, self.type]


@dataclass
class PhaseEventDef:
    """Holds code that defines some type of configuration of bodies that occurs once each period of a bound system."""

    name: str
    calculation: Callable[[str, Ephemeris, int], PhaseEvent]


def calc_time_of_phase(name: str, ephemeris: Ephemeris, orbit: int, phase: float) -> PhaseEvent:
    """Calculates events of a bound system that occur at a particular point in the phase of the system

    :param name: Descriptive title of the event, such as "mid eclipse"
    :type name: str
    :param ephemeris: Parameters of the bound system
    :type ephemeris: Ephemeris
    :param orbit: Orbit number in which the calculation should occur, starting with orbit zero at t0
    :type orbit: int
    :param phase: The desired phase value, from 0. to 1., at which the event occurs
    :type phase: float"""
    return PhaseEvent(
        time=Time(ephemeris.t0 + orbit * ephemeris.period + phase * ephemeris.period, format="jd"),
        system=ephemeris.system,
        member=ephemeris.component,
        orbit=orbit,
        type=name,
        phase=phase,
    )


def calc_mid_eclipse(name: str, ephemeris: Ephemeris, orbit: int) -> PhaseEvent:
    return calc_time_of_phase(name, ephemeris, orbit, 0.0)


def calc_time_of_gress(
    name: str, ephemeris: Ephemeris, orbit: int, ingress: bool = True, synthetic_phase_percent: float = float("NaN")
) -> PhaseEvent:
    """Calculates the time of ingress or egress of an eclipsing system

    :param name: Descriptive title of the event, such as "egress"
    :type name: str
    :param ephemeris: Parameters of the bound system
    :type ephemeris: Ephemeris
    :param orbit: The orbit number in which the calculation should occur, starting with orbit zero at t0
    :type orbit: int
    :param ingress: If true, calculate time of ingress into eclipse, otherwise, calculate time of egress
    :type ingress: bool
    :param synthetic_phase_percent: If duration is missing from ephemeris, pretend eclipse lasts this percentage of the period
    :type synthetic_phase_percent: float
    """
    t = ephemeris.t0 + orbit * ephemeris.period  # mid eclipse
    if ephemeris.duration == ephemeris.duration:
        half_duration = ephemeris.duration / 2
    else:
        if synthetic_phase_percent == synthetic_phase_percent:
            half_duration = ephemeris.period * synthetic_phase_percent / 2
        else:
            half_duration = u.Quantity(0, unit=u.hour)
    phase_change = (half_duration / ephemeris.period).to(u.dimensionless_unscaled)
    if ingress:
        t += ephemeris.period - half_duration
        phase = 1 - phase_change
    else:
        t += half_duration
        phase = phase_change
    return PhaseEvent(time=t, system=ephemeris.system, member=ephemeris.component, orbit=orbit, type=name, phase=phase)


class PhaseEventList:
    """Holds a list of PhaseEvent instances that occur between a specific start and end time"""

    def __init__(
        self,
        beg: Time = None,
        end: Time = None,
        beg_phase: float = float("nan"),
        end_phase: float = float("nan"),
        events: list[PhaseEvent] = None,
    ):
        self.beg = beg
        self.end = end
        self.beg_phase = beg_phase
        self.end_phase = end_phase
        if events:
            self.events = events
        else:
            self.events = []

    def __repr__(self):
        answer = f"{self.beg} to {self.end} \n"
        answer += f"Phase at start: {self.beg_phase:.2f}\n"
        answer += f"Phase at end: {self.end_phase:.2f}\n"
        if self.events and len(self.events) > 0:
            for event in self.events:
                answer += f"{event}\n"
        else:
            answer += "No events\n"
        return answer

    def calc_time_in_spans(self, beg: float, end: float) -> defaultdict[str, float]:
        """Calculate which event occupied the largest portion of a subsegment of an observing session"""
        if beg < self.beg or end > self.end or end < beg:
            raise ValueError(f"Parameters must be ordered, and within the timeframe {self.beg} to {self.end}")
        if not self.events:
            raise ValueError("This list has no events")
        time_spans = defaultdict(float)
        prev_event = None
        for event in self.events:
            if not event.time or event.tim < beg:
                # initial starting state is given by either most recent event before `beg` or a nan event
                prev_event = replace(event, time=beg)
                continue
            if event.time >= end:
                break  # there might be events after this time sequence, but we don't need them here
            if beg <= event.time < end:
                if prev_event:
                    this_span = min(event.time, end) - prev_event.time
                    time_spans[prev_event.type] = this_span
            prev_event = event
        if prev_event and prev_event.time < end:
            # capture the time we missed when the loop terminated, ie, from last event to `end`
            time_spans[prev_event.type] += end - prev_event.time
        if prev_event and len(time_spans) == 0:
            # if no events other than the one indicating starting state, we were in starting state the whole time
            time_spans[prev_event.type] = end - beg
        return time_spans

    def calc_longest_span(self, beg: float, end: float) -> str:
        spans = self.calc_time_in_spans(beg, end)
        return max(spans, key=spans.get)

    @staticmethod
    def calc_phase_events(ephem: Ephemeris, event_defs: list[PhaseEventDef], beg: Time, end: Time) -> "PhaseEventList":
        """Find all phase events `beg` <= event < `end`
        If no event exactly at `beg`, first returned event is last event before `beg`, with `time` None & `phase` set to `nan`
        The `event_defs` elements MUST be arranged in order of increasing phase
        """

        def calc_phase(t0: Time, period: u.Quantity, time: Time) -> float:
            orbit = int((time - t0) / period.to(u.day))
            last_eclipse = t0 + period * orbit
            return ((time - last_eclipse) / period).to(u.dimensionless_unscaled)

        if end < beg:
            raise ValueError("Time segment not chronological")
        if not ephem or not event_defs:
            raise ValueError("No parameters can be None")
        orbit = int((beg - ephem.t0) / ephem.period) - 1  # previous orbit (the -1) helps when window is short & event defs are few
        prev_event = None
        i = 0
        events = []
        while True:
            event_def = event_defs[i]
            event = event_def.calculation(event_def.name, ephem, orbit)
            if event.time < beg:
                # we only want to keep the most recent event before beg, to signify the state when the time region began
                # nullify the time and phase so we only return events beg <= event < end
                events = [replace(event, time=None, phase=float("nan"))]
            elif beg <= event.time < end:
                # these are the events we will keep
                if event.time == beg:
                    # clear out any previously kept starting event since we now have an event exactly at the start
                    events = []
                if not prev_event or prev_event.type != event.type:
                    # keep the event if it is of different type than the previous one
                    events.append(event)
            elif event.time >= end:
                break
            i = (i + 1) % len(event_defs)
            if i == 0:
                orbit += 1
            prev_event = event
        if len(events) == 0 and prev_event:
            # if we only had an event before beg and none between beg & end, record that prior event as in effect the whole time
            events = [replace(prev_event, time=None, phase=float("nan"))]
        beg_phase = calc_phase(ephem.t0, ephem.period, beg)
        end_phase = calc_phase(ephem.t0, ephem.period, end)
        answer = PhaseEventList(beg=beg, end=end, beg_phase=beg_phase, end_phase=end_phase, events=events)
        return answer
