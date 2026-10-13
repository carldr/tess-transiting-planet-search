"""Write a results CSV one row at a time, so that a stopped run can continue.

IncrementalCsvWriter writes each row to disk as soon as the step has it. A step reads the rows of its earlier run
with earlier_rows_to_keep, gives them to the writer, and skips their items. --restart makes the step start again.
"""

import csv
import os


def complete_csv_rows(path):
    """Return the rows of an earlier CSV. A stopped run can leave a partial last line, so the function drops it."""
    if not path.exists():
        return []
    text = path.read_text()
    if not text.endswith("\n"):
        text = text[:text.rfind("\n") + 1]
    return list(csv.DictReader(text.splitlines()))


def earlier_rows_to_keep(path, restart):
    """Return no rows with --restart. Otherwise, return the complete rows of the earlier run that have no error.

    The step does the items of the other rows again, so it retries a temporary failure, such as an archive outage.
    """
    if restart:
        return []
    return [row for row in complete_csv_rows(path) if not row.get("error")]


class IncrementalCsvWriter:
    """Write the rows of an earlier run, then each new row, and sync the file to disk after each row.

    The writer first writes the earlier rows to a temporary file and moves it over the CSV. A stop at that time
    thus cannot lose the earlier rows.
    """

    def __init__(self, path, columns, earlier_rows=()):
        temporary_path = path.with_name(path.name + ".partial")
        with open(temporary_path, "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(earlier_rows)
            file.flush()
            os.fsync(file.fileno())
        os.replace(temporary_path, path)
        self.file = open(path, "a", newline="")
        self.writer = csv.DictWriter(self.file, fieldnames=columns, extrasaction="ignore")

    def write_row(self, row):
        """Write the row, and return it as a later run reads it back: the text of each column, "" for None."""
        written_row = {column: "" if row.get(column) is None else str(row.get(column, ""))
                       for column in self.writer.fieldnames}
        self.writer.writerow(written_row)
        self.file.flush()
        os.fsync(self.file.fileno())
        return written_row

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.file.close()
