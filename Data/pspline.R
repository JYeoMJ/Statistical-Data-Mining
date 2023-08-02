pspline<-function(x, y, xnew=x, knots=5)
{
n = length(y)

orderx = order(x)
J = x[orderx]
K = floor(n/knots);
K = (1:(knots-1))*K
knotsp = J[K];
J = length(knotsp)
X = matrix(1, n, J+4)
X[,2] = x
X[,3] = x^2
X[,4] = x^3

ne = length(xnew)
Xe = matrix(1, ne, J+4)
Xe[,2] = xnew
Xe[,3] = xnew^2
Xe[,4] = xnew^3

for (j in 1:J)
{
   X[,4+j] = (x-knotsp[j])^3*(x-knotsp[j]>0)
   Xe[,4+j] = (xnew-knotsp[j])^3*(xnew-knotsp[j]>0)
}

invXX = solve(t(X) %*% X + diag(c(matrix(1, J+4, 1)))/n^2)
theta = invXX %*% (t(X) %*% y)

s2 = mean((y - X%*% theta)^2)

W = sqrt(diag(Xe %*% (invXX*s2) %*% t(Xe)))

m = Xe %*% theta


return(list(m=m, Un = m+1.96*W, Ln = m-1.96*W))
}