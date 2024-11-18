from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass, replace
from functools import partial


@dataclass
class ephemeris:
    system: str
    component: str
    t0: float
    period: float
    duration: float
    eccentricity: float


@dataclass
class phaseEvent:
    jd: float # a value of nan means the event occurred sometime before the current time segment
    system: str
    component: str
    orbit: int
    type: str
    phase: float

    def copy(self) -> "phaseEvent":
        return phaseEvent(self.jd, self.system, self.component, self.orbit, self.type, self.phase)

# helper functions useful when making phaseEventDef instances
def calc_time_of_phase(name: str, ephemeris: ephemeris, orbit: int, phase: float) -> list[phaseEvent]:
    return phaseEvent(
        jd=ephemeris.t0 + orbit * ephemeris.period + phase * ephemeris.period,
        system=ephemeris.system,
        component=ephemeris.component,
        orbit=orbit,
        type=name,
        phase=phase,
    )


def calc_mid_eclipse(name: str, ephemeris: ephemeris, orbit: int) -> list[phaseEvent]:
    return calc_time_of_phase(name, ephemeris, orbit, 0.0)


def calc_time_of_gress(name: str, ephemeris: ephemeris, orbit: int, ingress: bool = True) -> list[phaseEvent]:
    t = ephemeris.t0 + orbit * ephemeris.period  # mid eclipse
    half_duration = ephemeris.duration / 2
    phase_change = half_duration / ephemeris.period
    if ingress:
        t += ephemeris.period - half_duration
        phase = 1 - phase_change
    else:
        t += half_duration
        phase = phase_change
    return phaseEvent(jd=t, system=ephemeris.system, component=ephemeris.component, orbit=orbit, type=name, phase=phase)


@dataclass
class phaseEventDef:
    name: str
    calculation: Callable[[ephemeris, int], phaseEvent]

class phaseEventList:

    def __init__(
        self,
        beg: float = float("nan"),
        end: float = float("nan"),
        beg_phase: float = float("nan"),
        end_phase: float = float("nan"),
        events: list[phaseEvent] = [],
    ):
        self.beg = beg
        self.end = end
        self.beg_phase = beg_phase
        self.end_phase = end_phase
        self.events = events

    def __repr__(self):
        answer = f"{self.beg} JD to {self.end} JD\n"
        answer += f"Phase at start: {self.beg_phase:.2f}\n"
        answer += f"Phase at end: {self.end_phase:.2f}\n"
        for event in self.events:
            answer += f"{event}\n"
        return answer

    def determine_phase_for_span(self, beg: float, end: float):
        """Calculate which event occupied the largest portion of a subsection of an observing session"""
        if beg < self.beg or end > self.end:
            raise ValueError(f"Parameters need to fall within the timeframe {self.beg} to {self.end}")
        time_spans = defaultdict(float)
        for event in self.events:
            if event.jd != event.jd:
                # first event is often jd=nan to indicate it occurred before start of this time sequence
                # treat it here as if it occurred at the beginning of this time sequence
                prev_event = replace(event, jd=beg)
                continue
            this_span = min(event.jd, end) - prev_event.jd
            time_spans[prev_event.type] = this_span
            if event.jd >= end:
                break # there might be events after this time sequence, but we don't need them here
            prev_event = event
        last_event = self.events[-1]
        if last_event.jd < end:
            #capture the time we missed when the loop terminated
            time_spans[last_event.type] += end - last_event.jd
        if len(time_spans) == 0:
            time_spans[prev_event.type] = end - beg
        return time_spans
        # return max(time_spans, key=time_spans.get)

    @staticmethod
    def calc_phase_events(ephem: ephemeris, eventDefs: list[phaseEventDef], beg: float, end: float) -> "phaseEventList":
        orbit = int((beg - ephem.t0) / ephem.period)
        i = 0
        t = ephem.t0 + orbit * ephem.period  # start at the mid-eclipse prior to beg
        if t >= beg:
            orbit -= 1
        event_at_beg = None
        def calc_phase(t0:float, period:float, time:float) -> float:
            orbit = int((time - t0) / period)
            last_eclipse = t0 + period * orbit
            return (time - last_eclipse) / period
        beg_phase = calc_phase(ephem.t0, ephem.period, beg)
        end_phase = calc_phase(ephem.t0, ephem.period, end)
        answer = phaseEventList(beg=beg, end=end, beg_phase=beg_phase, end_phase=end_phase)
        while t < end:
            span = eventDefs[i]
            event = span.calculation(span.name, ephem, orbit)
            if beg <= event.jd <= end:
                if len(answer.events) == 0:
                    # insert a possibly partial event to note the event in force at start of time interval
                    first_event = event_at_beg or event
                    if first_event.jd < beg:
                        first_event.jd = float("nan")
                        answer.events.append(first_event)
                    answer.events.append(event)
                else:
                    answer.events.append(event)
            else:
                event_at_beg = event
            i += 1
            if i == len(eventDefs):
                i = 0
                orbit += 1
            t = event.jd
        return answer