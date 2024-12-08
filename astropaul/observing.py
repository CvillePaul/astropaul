import operator

import astroplan as ap
from astropy.time import Time
import astropy.units as u
import pandas as pd


class ObservingSession:
    def __init__(
        self,
        observer: ap.Observer,
        observing_segments: list[tuple[Time, Time]] = None,
    ):
        self.observer = observer
        self.observing_segments = observing_segments or []

    def __repr__(self):
        if self.observer.name:
            answer = f"{self.observer.name} "
        else:
            answer = "Site "
        answer += f"({self.observer.latitude.value:+.3f}, {self.observer.longitude.value:+.3f})"
        answer += f" in {self.observer.timezone}\n"
        if self.observing_segments:
            answer += "Observing segments:\n"
            for beg, end in self.observing_segments:
                answer += f"    {beg.iso[:19]} to {end.iso[:19]}\n"
        else:
            answer += "No observing segments defined\n"
        return answer

    @property
    def time_range(self) -> tuple[Time, Time]:
        return (
            min(self.observing_segments, key=operator.itemgetter(0))[0],
            max(self.observing_segments, key=operator.itemgetter(1))[1],
        )

    def _determine_nighttime(self, night: Time) -> tuple[Time, Time]:
        beg_night = self.observer.sun_set_time(night, which="nearest")
        end_night = self.observer.sun_rise_time(beg_night, which="next")
        return (beg_night, end_night)

    def add_full_day(self, day: str | Time):
        self.observing_segments.append(self._determine_nighttime(Time(day)))

    def add_half_day(self, day: str | Time, first_half: bool = True):
        beg, end = self._determine_nighttime(Time(day))
        mid = Time(beg.jd + ((end.jd - beg.jd) / 2), format="jd")
        if first_half:
            self.observing_segments.append((beg, mid))
        else:
            self.observing_segments.append((mid, end))

    def add_day_range(self, beg: str | Time, end: str | Time):
        beg, end = Time(beg), Time(end)
        day = beg
        while day < end:
            self.add_full_day(day)
            day += 1 * u.day

    def calc_subsegments(
        self, interval: u.Quantity = 1 * u.hour, top_of_hour: bool = True, skip_partial: bool = True
    ) -> list[list[(Time, Time)]]:
        """Chop each observing segment into a list of subsegments that are `interval` long.
        In general, the last sub segment will be shorter than `interval`
        If `top_of_hour`, first sub segment might also be shorter than `interval`
        Sub segments shorter than `interval` can be skipped by setting `skip_partial`"""

        answer = []
        for beg, end in self.observing_segments:
            sub_segment = []
            t = beg
            if top_of_hour:
                t = Time(pd.Timestamp(beg.to_datetime()).ceil("h"))
                if t - beg > 0 * u.day and not skip_partial:
                    sub_segment.append((beg, t))  # 1st sub segment, smaller than `interval`
            while t < end:
                next = t + interval
                if next < end:
                    sub_segment.append((t, next))
                else:
                    if not skip_partial:
                        sub_segment.append((t, end))
                t = next
            answer.append(sub_segment)
        return answer
