import numpy as np

##############################################################################

def transmat(n,p,q):
    if (n <= 1):
        sys.exit('n needs to be bigger than 1 (transmat)')
    P0 = np.array([[p,1-p],[1-q,q]])
    for i in np.arange(2,n):
        Pnew = np.zeros([i+1,i+1])
        z0 = np.zeros([i+1,i+1])
        z1 = np.zeros([i+1,i+1])
        z2 = np.zeros([i+1,i+1])
        z3 = np.zeros([i+1,i+1])
        z0[:i,:i] = P0
        z1[:i,1:] = P0
        z2[1:,:i] = P0
        z3[1:,1:] = P0
        mat_sum = p*z0+(1-p)*z1+(1-q)*z2+q*z3
        Pnew[0,:] = mat_sum[0,:]
        Pnew[i,:] = mat_sum[i,:]
        Pnew[1:i,:] = 0.5*mat_sum[1:i,:]
        P0 = Pnew
    return(P0)

def solve_case(eta,P,params,ns,js):

    beta = params[0]
    gamma_pi = params[1]
    gamma_y = params[2]
    kappa = params[3]
    dpbar = params[4]
    Rss = np.log(dpbar/beta)

    eyens = np.identity(ns)
    H0 = np.zeros([2*ns,2*ns])
    H1 = np.ones([ns,ns])
    eyejs = np.identity(js)
    eyejps = np.identity(ns-js)

    H0 = np.zeros([2*ns,2*ns])
    H0[:js,:js] = eyejs
    H0[js:ns,js:ns] = (1+gamma_y)*eyejps
    H0[js:ns,ns+js:] = gamma_pi*eyejps
    H0[ns:,:ns] = -kappa*eyens
    H0[ns:,ns:] = eyens
    H1 = np.zeros([2*ns,2*ns])
    H1[:ns,:ns] = P
    H1[:ns,ns:] = P
    H1[ns:,ns:] = beta*P

    A = H0-H1
    B = np.zeros(2*ns)
    B[:js] = Rss+eta[:js]
    B[js:ns] = eta[js:ns]
    xx = np.linalg.solve(A,B)
    yy = xx[:ns]
    dp = xx[ns:]
    notr = gamma_pi*dp+gamma_y*yy
    nomr = np.maximum(-Rss,notr)
    return(yy,dp,notr,nomr)           

def solve_model(eta,P,params,ns):

    Rss = np.log(params[4]/params[0])
    for j in np.arange(ns+1):
        (yy,dp,notr,nomr) = solve_case(eta,P,params,ns,j)
        if (all(notr[:j]<=-Rss) and all(notr[j:] >= -Rss)):
            solution_type = j
            break
        else:
            solution_type = -1
    return(solution_type,yy,dp,notr,nomr)    
        
















    





