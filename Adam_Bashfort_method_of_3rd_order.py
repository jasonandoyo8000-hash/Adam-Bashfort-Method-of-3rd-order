#Adam-Bashfort method of 3rd order for numerical solutions of ordinary differential equations
from numpy import *
from matplotlib.pyplot import *



def AB_3rd(f,x0,xf,yn,h,r):
    """
    This function calculates the approximate values of the solution of the ODE with initial value y(x0) = yn
    in the interval [x0,xf] by using Adams-Bashforth Method - 3rd Order
    
    Input: AB_3rd(dy/dx, x0, xf, yn, h, r)
    where
    dy/dx is the ODE (usually in the form dy/dx = y'(x) = f(x,y))
    y(x0) = yn is the initial value
    h is the step size
    r is the number of decimals in the approximation
    
    Output: Two lists: the first list is compose of the approximate values of y and
    the second list is compose of the values of x
    """
    
    def arccos(x):
        return acos(x)


    x_list = [x0]
    first_3_values = [yn]
    
    #Runge-Kutta Method
    for xn in arange(x0,x0+2*h,h):
        k1 = f(xn,yn)
        k2 = f(xn+0.5*h,yn+0.5*h*k1)
        k3 = f(xn+0.5*h,yn+0.5*h*k2)
        k4 = f(xn+h,yn+h*k3)
        y = yn+(h/6)*(k1+2*k2+2*k3+k4)
        x_list.append(round(xn+h,r))
        first_3_values.append(round(y,r))
        yn = y

    #Adam-Bashfort Method - 3rd Order
    list1 = [first_3_values[0], first_3_values[1], first_3_values[2]]
    y1, y2, y3 = first_3_values[0], first_3_values[1], first_3_values[2]

    for xn in arange(x0,xf-2*h,h):
        y_ab = h*((5/12)*f(xn,y1)-(4/3)*f(xn+h,y2)+(23/12)*f(xn+2*h,y3))+y3
        list1.append(round(y_ab, r))
        x_list.append(round(xn+3*h,r))
        y1 = y2
        y2 = y3
        y3 = y_ab
        
    return list1, x_list

def disp(list1,list2):
    """
    This function displays the approximate values of the given ODE calculated by the function AB_3rd
    in matrix form where the first column are the values of x and the second column are the approximate
    values of y
    """
    
    print("First column are the values of x and the second column are the approximated values of y.")
    list3 = []
    
    for i in range(0,len(list1)):
        list3.append([list1[i],list2[i]])
    return array(list3)


def slope(f, list1, list2):
    """
    This function graph the approximated value of the ODE calculated by the function AB_3rd
    and shows the direction field of the given ODE
    """
    def arccos(X):
        list1=[]
        for j in range(0,len(X)):
            list2 = []
            for i in X[j]:
                list2.append(acos(i))
            list1.append(list2)
        return array(list1)

    title("The graph of approximated values of the solution of the given ODE")
    x = arange(min(list1),max(list1),(max(list1)-min(list1))/20)
    y = arange(min(list2),max(list2),(max(list2)-min(list2))/20)

    X,Y = meshgrid(x,y)

    dy = f(X,Y)
    dx = ones(dy.shape)

    norm = sqrt(dy**2+dx**2)

    dyu = dy/norm

    dxu = dx/norm
    
    quiver(X,Y,dxu,dyu,color = "r")
    xlim(min(list1),max(list1))
    ylim(min(list2),max(list2))
    plot(list1,list2, "k")
    grid()
    show()

def drive(a,b,c,d,e,p):
    """
    This function is for summarizing all the functions above
    """
    g = AB_3rd(a,b,c,d,e,p)

    k = slope(a,g[1],g[0])

    h = disp(g[1],g[0])
    
    return h



print("""Assuming that we have an ODE (Ordinary Differential Equation) in the form

                              dy/dx = f(x,y)
                              
with initial value y(x0)=y0. This python code will give the approximate values of 
the solution of the above IVP (Initial Value Problem) on the interval [x0, xf]
using Adams-Bashforth Method - 3rd Order with step size h. In addition, the direction
field of the ODE will also be presented along with the graph of the approximated values
of the solution of the IVP. We denote r as the number of decimals. 


IMPORTANT NOTE: In typing for f(x,y), please use the syntax of python, for example, 
instead of typing x^2, please type x**2. 


Enter the following: """)


dydx = input("f(x,y) = ")
x0 = eval(input("x0 = "))
xf = eval(input("xf = "))
y0 = eval(input("y0 = "))
h = eval(input("h = "))
r = eval(input("r = "))

exec('f = lambda x, y :'+dydx) 

print(drive(f, x0, xf, y0, h, r))
