"""Write a results CSV one row at a time, so that a stopped run can continue.

IncrementalCsvWriter writes each row to disk as soon as the step has it. A step that runs with --resume reads the
rows of its earlier run with complete_csv_rows, gives them to the writer, and skips them.
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
        self.writer.writerow(row)
        self.file.flush()
        os.fsync(self.file.fileno())

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.file.close()
