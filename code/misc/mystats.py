import math
import numpy.random


def mean(values):
    values = filter(lambda x: x != 'NA' and not math.isnan(x), values)
    if len(values) == 0:
        # TODO float('nan') ?
        return 'NA'
    return float(sum(values)) / len(values)


def std_dev(values):
    values = filter(lambda x: x != 'NA' and not math.isnan(x), values)
    if len(values) == 0:
        return 'NA'
    if len(values) == 1:
        return 0
    m = mean(values)
    return math.sqrt(sum([(x - m)**2 for x in values]) / (len(values) - 1))


def std_err(values):
    values = filter(lambda x: x != 'NA' and not math.isnan(x), values)
    if len(values) == 0:
        return 'NA'
    return std_dev(values) / math.sqrt(len(values))


def bootstrap(values, n=100, alpha=.05):
    values = filter(lambda x: x != 'NA' and not math.isnan(x), values)
    x = len(values)
    if x == 0:
        return 'NA', 'NA'
    a = []
    for i in range(n):
        a.append(mean(numpy.random.choice(values, size=x, replace=True)))
    a.sort()
    # print len(a), a.count(0)
    # print mean(a)
    return a[int(alpha * n * .5)], a[int((1 - alpha * .5) * n)]


def median(values):
    m = sorted(values)
    x = len(m)
    if x % 2 == 0:
        return mean([m[x/2], m[x/2-1]])
    return m[(x-1)/2]
