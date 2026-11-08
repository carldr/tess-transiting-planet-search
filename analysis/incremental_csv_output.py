"""Write a results CSV one row at a time, so that a stopped run can continue.

IncrementalCsvWriter writes each row to disk as soon as the step has it. A step reads the rows of its earlier run
with earlier_rows_to_keep, gives them to the writer, and skips their items. --restart makes the step start again.
--tic-id limits a step to the candidates of the given TIC IDs, and limits --restart to their rows.
With compress_when_finished, a finished run leaves only a gzip copy of the CSV, and the next run reads its rows.
"""

import csv
import gzip
import os
import shutil


def compressed_path(path):
    """Return the path of the gzip copy that IncrementalCsvWriter leaves of path when compress_when_finished is set."""
    return path.with_name(path.name + ".gz")


def read_csv_rows(path):
    """Return the rows of a CSV, read through gzip when the name ends in .gz."""
    with (gzip.open(path, "rt", newline="") if path.suffix == ".gz" else open(path, newline="")) as file:
        return list(csv.DictReader(file))


def complete_csv_rows(path):
    """Return the rows of an earlier CSV. A stopped run can leave a partial last line, so the function drops it.

    When path does not exist, the function reads the gzip copy of a finished run, if there is one.
    """
    if not path.exists():
        return read_csv_rows(compressed_path(path)) if compressed_path(path).exists() else []
    text = path.read_text()
    if not text.endswith("\n"):
        text = text[:text.rfind("\n") + 1]
    return list(csv.DictReader(text.splitlines()))


def earlier_rows_to_keep(path, restart, input_keys, key_columns=("tic_id", "period_days"), tic_ids=()):
    """Return the complete rows of the earlier run that have no error and whose item is still in the input.

    input_keys holds the key_columns values of each item in the current input of the step. An earlier row of an item
    that is no longer in the input, such as a candidate that now matches a TOI, has no place in the output.
    The step does the items of the other rows again, so it retries a temporary failure, such as an archive outage.
    With --restart, the function also drops the rows of tic_ids, or every row when tic_ids is empty.
    """
    rows = [row for row in complete_csv_rows(path)
            if not row.get("error") and tuple(row[column] for column in key_columns) in input_keys]
    if restart:
        return [row for row in rows if tic_ids and row["tic_id"] not in tic_ids]
    return rows


class IncrementalCsvWriter:
    """Write the rows of an earlier run, then each new row, and sync the file to disk after each row.

    The writer first writes the earlier rows to a temporary file and moves it over the CSV. A stop at that time
    thus cannot lose the earlier rows.
    """

    def __init__(self, path, columns, earlier_rows=(), compress_when_finished=False):
        self.path, self.compress_when_finished = path, compress_when_finished
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

    def __exit__(self, exception_type, exception, traceback):
        self.file.close()
        # A finished run replaces the CSV with its gzip copy. A stopped run keeps the CSV, so the next run continues.
        if self.compress_when_finished and exception_type is None:
            temporary_path = compressed_path(self.path).with_name(compressed_path(self.path).name + ".partial")
            with open(self.path, "rb") as source, gzip.open(temporary_path, "wb") as target:
                shutil.copyfileobj(source, target)
            os.replace(temporary_path, compressed_path(self.path))
            self.path.unlink()
