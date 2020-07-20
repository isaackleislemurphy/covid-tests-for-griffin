
import numpy as np
import matplotlib.pyplot as plt

def tests_per_pool(p, n_max = 25):
    
    '''
    Calculates the number of tests per person for a given prevalence and various pool sizes. 
    Used for optimizing pool size. 
    
    @param <float> p: prevalence, \in [0, 1]
    @param <int> n_max: maximum pool size
    @return: a list of tests per person corresponding to each pool size for the prevalence
    '''
    
    phi = 1 - p
    
    tests_per_pool = [phi**i  + (i + 1)*(1 - phi**i) for i in range(2, n_max + 1)]
    tests_per_person =[(phi**i  + (i + 1)*(1 - phi**i))/i for i in range(2, n_max + 1)]
    
    plt.plot(range(2, n_max + 1), tests_per_person)
    plt.xlabel('N')
    plt.ylabel('Tests per Person')
    plt.ylim(0, 1.5)
    plt.plot(range(2, n_max + 1), [1 for i in range(2, n_max + 1)])
    plt.show()
    
    return tests_per_person



def newton(f, f_prime, x0 = .85, epsilon = 1e-6, max_iter = 1000, kwargs = [4]):
    
    '''
    A Newtons method root finder. 
    
    @param <func> f: the differentiable function to be rooted
    @param <func> f_prime: the derivative of f
    @param <float> x0: initial guess
    @param <float> epsilon: stopping tolerance
    @param <int> max_iter: maximum number of iterations to be performed
    @param <list> kwargs: items to be passed to f and/or f_prime
    
    @return: the root of the function, a float
    '''
    
    xn = x0
    for n in range(0,max_iter):
        fxn = f(xn, *kwargs)
        if abs(fxn) < epsilon:
            return xn
        dfxn = f_prime(xn, *kwargs)
        if dfxn == 0:
            print('Zero derivative. No solution found.')
            return None
        xn = xn - fxn/dfxn
    print('Exceeded maximum iterations. No solution found.')
    return None



def exp_tests(phi, n, breakeven = 1):
    
    '''
    Computes the number of expected tests per person, given a pool size, less a breakeven rate
    
    @param <float> phi: 1 - prevalence
    @param <int> n: pool size
    @param <float> breakeven: a tests per person breakeven rate, between 0 and 1
    
    @return: float, E[Y]/n - breakeven
    '''
    
    return (phi**n + (n + 1)*(1 - phi**n))/n - breakeven



def exp_tests_prime(phi, n, *kwargs):
    
    '''
    A function that computes the d/dPhi E[Y]. See exp_tests()
    
    @return: float, d/dPhi (E[Y]/n - breakeven)
    '''
    
    return (-n)*phi**(n - 1)



def solve_prevalence(n = 4, breakeven = 1, guess = .85):
    
    '''
    Main function. Given a n, and a breakeven rate, it identifies the prevalence necessary to achieve the
    breakeven rate. 
    
    @param <int> n: pool size
    @param <float> breakeven: the desired breakeven rate
    @param <float> guess: an initial phi guess for the Newton solver
    
    @return: a float, the prevalence that would acheive the desired breakeven. 
    
    '''
    
    return 1 - newton(exp_tests, exp_tests_prime, x0 = guess, kwargs = [n, breakeven])

    

def plot_breakeven(n = 4):
    
    rng = np.arange(0, 1, .001)
    plt.plot(rng, [solve_prevalence(n, breakeven  = i) for i in rng])    
    plt.ylim(0, )
    plt.xlabel('Desired Tests Per Person')
    plt.ylabel('Necessary Prevalence for Des. Tests/Person')
    plt.show()
