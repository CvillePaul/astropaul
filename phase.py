# phase event: jd type phase
# ephemeris: t0 period duration eccentricity
# phase span: begphase endphase name

# phase event calculator: ephemeris eventfuncs gresses phasespans
#   calc_phase_events: beg end
#   determine_phase_span: (phaseevents) beg end

from dataclasses import dataclass


@dataclass
class phaseEvent:
    jd: float
    type: str
    phase: float


@dataclass
class ephemeris:
    t0: float
    period: float
    duration: float
    eccentricity: float


@dataclass
class phaseSpan:
    beg_phase: float
    end_phase: float
    name: str


class phaseEventCalculator:
    def __init__(
        self,
        ephem: ephemeris,
        event_funcs: list,
        calc_gresses: bool,
        spans: list[phaseSpan],
    ):
        self.ephem = ephem
        self.event_funcs = event_funcs
        self.calc_gresses = calc_gresses
        self.spans = spans

    def calc_phase_events(self, beg: float, end: float) -> list[phaseEvent]:
        pass

    def determine_phase_span(self, beg: float, end: float) -> phaseSpan:
        pass
