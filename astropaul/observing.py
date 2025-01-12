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
        answer = self.site_info + "\n"
        if self.observing_segments:
            answer += f"{self.total_time.to(u.hour):.1f} total observing time\n"
            num_segments = len(self.observing_segments)
            if num_segments < 10:
                answer += f"Observing segments ({num_segments}):\n"
                for beg, end in self.observing_segments:
                    answer += f"    {(end - beg).to(u.hour).value:.1f} hours: {beg.iso[:19]} to {end.iso[:19]}\n"
            else:
                answer += f"{num_segments} segments, from {self.observing_segments[0][0].iso[:19]} to {self.observing_segments[-1][-1].iso[:19]}\n"
        else:
            answer += "No observing segments defined\n"
        return answer

    @property
    def site_info(self) -> str:
        if self.observer.name:
            answer = f"{self.observer.name} "
        else:
            answer = "Site "
        answer += f"({self.observer.latitude.value:+.3f}, {self.observer.longitude.value:+.3f})"
        answer += f" timezone {self.observer.timezone}"
        return answer

    @property
    def total_time(self) -> Time:
        answer = 0 * u.hour
        for beg, end in self.observing_segments:
            answer += end - beg
        return answer

    @property
    def time_range(self) -> tuple[Time, Time]:
        return (
            min(self.observing_segments, key=operator.itemgetter(0))[0],
            max(self.observing_segments, key=operator.itemgetter(1))[1],
        )

    @property
    def starting_lst(self) -> float:
        if not self.observing_segments:
            return 0
        else:
            return self.observer.local_sidereal_time(self.observing_segments[0][0]).to(u.deg).value

    def _determine_nighttime(self, night: Time) -> tuple[Time, Time]:
        beg_night = self.observer.sun_set_time(night, which="nearest")
        end_night = self.observer.sun_rise_time(beg_night, which="next")
        return (beg_night, end_night)

    def add_full_day(self, day: str | Time):
        self.observing_segments.append(self._determine_nighttime(Time(day)))
        return self

    def add_half_day(self, day: str | Time, first_half: bool = True):
        beg, end = self._determine_nighttime(Time(day))
        mid = Time(beg.jd + ((end.jd - beg.jd) / 2), format="jd")
        if first_half:
            self.observing_segments.append((beg, mid))
        else:
            self.observing_segments.append((mid, end))
        return self

    def add_day_range(self, beg: str | Time, end: str | Time):
        end_time = Time(end) + 1 * u.day
        day = Time(beg)
        while day < end_time:
            self.add_full_day(day)
            day += 1 * u.day
        return self

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
