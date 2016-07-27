# a utility function
my.f <- function(x) { sum(x[x>0]*log(x[x>0])) - sum(x)*log(sum(x)) }

# the test statistic D (as seen in the wiki article)
my.D <- function(x, y) { -2*(my.f(x+y) - my.f(x) - my.f(y)) }

# the pvalue. x and y are counts.
my.pval <- function(x, y) { pchisq(my.D(x,y), length(x)-1, lower.tail=F) }