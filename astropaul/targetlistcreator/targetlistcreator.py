import collections
import datetime
from sqlite3 import Connection
from typing import Any

import numpy as np
import pandas as pd

from astropaul.database import db_style_to_string

import __main__


class TargetListCriterion:
    def __init__(self, criterion: str, children: list["TargetListCriterion"] = None):
        self.criterion = criterion
        self.children = children or list()

    def copy(self) -> "TargetListCriterion":
        answer = TargetListCriterion(criterion=self.criterion, children=self.children.copy())
        return answer

    def __repr__(self) -> str:
        answer = self.criterion
        if self.children and len(self.children) > 0:
            for child in self.children:
                for line in str(child).split("\n"):
                    answer += f"\n  {line}"
        return answer


class TargetListCriteria:
    def __init__(self, criteria: list[TargetListCriterion] = None):
        self.criteria = criteria or []

    def copy(self) -> "TargetListCriteria":
        answer = TargetListCriteria([criterion.copy() for criterion in self.criteria])
        return answer

    def add(self, criterion) -> None:
        if isinstance(criterion, str):
            self.criteria.append(TargetListCriterion(criterion))
        elif isinstance(criterion, TargetListCriterion):
            self.criteria.append(criterion)
        else:
            raise ValueError("Bad type for criterion")

    def __repr__(self) -> str:
        answer = ""
        for criterion in self.criteria:
            answer += f"{criterion}\n"
        return answer

    def __len__(self) -> int:
        return len(self.criteria)

    def __iter__(self):
        return iter(self.criteria)


class TargetList:
    def __init__(
        self,
        name: str = "Target List",
        target_list: pd.DataFrame = None,
        list_criteria: TargetListCriteria = None,
        column_groups: dict[str, tuple[list[str], list[str]]] = None,
        other_lists: dict[str, pd.DataFrame] = None,
    ):
        self.name = name
        self.target_list = target_list if target_list is not None else pd.DataFrame()
        self.list_criteria = list_criteria or TargetListCriteria()
        self.column_groups = column_groups or collections.defaultdict(list)
        if not other_lists:
            self.other_lists = dict()
        else:
            self.other_lists = other_lists

    def copy(self) -> "TargetList":
        answer = TargetList(
            name=self.name,
            target_list=self.target_list.copy(),
            list_criteria=self.list_criteria.copy(),
            column_groups=self.column_groups.copy(),
            other_lists={key: value.copy() for key, value in self.other_lists.items()},
        )
        return answer

    def add_columns(self, columns: dict[str, Any] = {}):
        pass

    def add_other(self, name: str, other: pd.DataFrame):
        self.other_lists[name] = other

    def summarize(self) -> str:
        answer = f"Name: {self.name}\n"
        answer += "Criteria:\n"
        if len(self.list_criteria) > 0:
            for line in str(self.list_criteria).split("\n"):
                answer += f"  {line}\n"
        else:
            answer += "    (none)\n"
        target_types = collections.Counter(self.target_list["Target Type"])
        answer += f"{len(self.target_list)} targets:\n"
        for type, count in target_types.items():
            answer += f"  {count:4d} {type}\n"
        answer += "Column Count (primary, secondary):\n"
        for name, (primary, secondary) in self.column_groups.items():
            answer += f"    {name}: ({len(primary)}, {len(secondary)})\n"
        answer += "Associated tables:\n"
        for name, other in self.other_lists.items():
            answer += f"    {len(other):4d} rows, {len(other.columns):2d} columns: {name}\n"
        if "PEPSI exp_time" in self.target_list.columns:
            answer += f"\nTotal PEPSI exposure time: {np.sum(self.target_list["PEPSI exp_time"])/60:.1f} minutes"
        return answer

    @staticmethod
    def union(tls: list["TargetList"], name: str = "Merged Target List") -> "TargetList":
        if not tls or len(tls) == 0:
            raise ValueError("Must supply at least one TargetList")
        answer = tls[0].copy()
        answer.name = name
        for tl in tls[1:]:
            answer.target_list = pd.concat([answer.target_list, tl.target_list], ignore_index=True)
            answer.column_groups |= tl.column_groups
            for table_name, table in tl.other_lists.items():
                if table_name in answer.other_lists:
                    answer.other_lists[table_name] = pd.concat([answer.other_lists[table_name], table])
                else:
                    answer.other_lists[table_name] = table
        answer.list_criteria = TargetListCriteria(
            [
                TargetListCriterion(
                    f"Union of lists: {" and ".join([f"{tl.name}" for tl in tls])}:",
                    TargetListCriteria([TargetListCriterion(f"List {tl.name}:", tl.list_criteria) for tl in tls]),
                )
            ]
        )
        return answer


class TargetListCreator:
    def __init__(self, name: str = "Standard", connection: Connection = None, steps: list = None, **kwargs):
        self.name = name
        self.connection = connection
        self.kwargs = kwargs.copy()
        self.steps = steps

    def calculate(
        self, initial_list: TargetList = None, steps=None, name: str = None, verbose: bool = False, **kwargs
    ) -> TargetList:
        """
        Create a new target list by running through self.steps and returning the result
        """
        intermediate_tl = initial_list.copy() if initial_list else TargetList(name=name or self.name)
        merged_kwargs = {"connection": self.connection, **self.kwargs, **kwargs}
        if steps is None:
            steps = self.steps
        for step in steps:
            intermediate_tl = step(intermediate_tl, **merged_kwargs)
            if verbose:
                print(f"{len(intermediate_tl.target_list):4d} targets, {datetime.datetime.now()}, {step=}")
        if verbose:
            print(intermediate_tl.summarize())
        return intermediate_tl


# some convenience functions
def convert_columns_to_human_style(df: pd.DataFrame) -> None:
    if index_name := df.index.name:
        df.index.name = db_style_to_string(index_name)
    df.columns = map(db_style_to_string, df)


def verify_step_requirements(tl: TargetList, required_tables: set[str] = None, required_columns: set[str] = None):
    if tl.target_list.empty:
        raise ValueError("Target List cannot be empty")
    if required_tables:
        existing_tables = set(tl.other_lists.keys())
        if not required_tables.issubset(existing_tables):
            raise ValueError(f"Required table(s) missing: {required_tables - existing_tables}")
    if required_columns:
        existing_columns = set(tl.target_list.columns)
        if not required_columns.issubset(existing_columns):
            raise ValueError(f"Required column(s) missing: {required_columns - existing_columns}")


# following are methods are meant to be used as steps of a TargetListCreator
