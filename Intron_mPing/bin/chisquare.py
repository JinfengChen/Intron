import numpy as np
from scipy.stats.mstats import chisquare
from scipy.stats import fisher_exact

observed = np.array([34, 111])
expected = np.array([71, 281])
a = chisquare(observed, expected)
b = fisher_exact([observed, expected])

print a
print b
