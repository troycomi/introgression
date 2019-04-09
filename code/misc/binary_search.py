import bisect
from typing import List


def present(a: List[int], x: int) -> bool:
    'Locate the leftmost value exactly equal to x in a'
    i = bisect.bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return True
    return False
