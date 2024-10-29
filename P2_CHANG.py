import argparse
from ortools.linear_solver import pywraplp
import numpy as np

def solve_with_ortools(constraints, objective_coeffs, maximize=True):   
    solver = pywraplp.Solver.CreateSolver("GLOP")
    if not solver:
        return

    x = solver.NumVar(0, solver.infinity(), "x")
    y = solver.NumVar(0, solver.infinity(), "y")

    print("Number of variables =", solver.NumVariables())

    for i, (coeff_x, coeff_y, bound) in enumerate(constraints):
        solver.Add(coeff_x * x + coeff_y * y <= bound)
    
    print("Number of constraints =", solver.NumConstraints())

    coeff_x, coeff_y = objective_coeffs
    if maximize:
        solver.Maximize(coeff_x * x + coeff_y * y)
    else:
        solver.Minimize(coeff_x * x + coeff_y * y)

    print(f"Solving with {solver.SolverVersion()}")
    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:
        print("Solution:")
        print(f"Objective value = {solver.Objective().Value():0.1f}")
        print(f"x = {x.solution_value():0.1f}")
        print(f"y = {y.solution_value():0.1f}")

        # Printing constraints
        print("\nVerification of constraints:")
        for i, (coeff_x, coeff_y,bound) in enumerate(constraints):
            constraint_value = coeff_x * x.solution_value() + coeff_y * y.solution_value()
            print(f"Constraint {i}: {coeff_x}*x + {coeff_y}*y <= {bound}")
    else:
        print("The problem does not have an optimal solution.")

    print("\nAdvanced usage:")
    print(f"Problem solved in {solver.wall_time():d} milliseconds")
    print(f"Problem solved in {solver.iterations():d} iterations")

def print_tableau(tableau):
    print("\nCurrent tableau :")
    num_rows, num_cols = tableau.shape
    header = ["Basic Var"] + ['x', 'y'] + [f's{i + 1}' for i in range(num_rows - 1)] + ['RHS']
    print("{:<15} {}".format(header[0], "  ".join(f"{h:<10}" for h in header[1:])))
    
    for i in range(num_rows):
        row = tableau[i]
        basic_var = f"x_{i + 1}" if i < num_rows - 1 else "Z"
        print(f"{basic_var:<15}", end=" ")
        print("  ".join(f"{val:<10.2f}" for val in row))
    print()

def simplex_step_by_step(c, A, b):
    num_vars = len(c)
    num_constraints = A.shape[0]

    tableau = np.zeros((num_constraints + 1, num_vars + num_constraints + 1))
    
    # Fill the tableau
    tableau[:-1, :num_vars] = A  # Coefficients of constraints
    tableau[:-1, num_vars:num_vars + num_constraints] = np.eye(num_constraints)  # Basic variables
    tableau[:-1, -1] = b  # RHS

    # Coefficients of the objective function
    tableau[-1, :num_vars] = -c

    print_tableau(tableau)

    while True:
        # Entring variable
        if np.all(tableau[-1, :-1] >= 0):
            print("Optimal solution achieved")
            break
        
        pivot_col = np.argmin(tableau[-1, :-1])  
        if np.all(tableau[:-1, pivot_col] <= 0):
            print("The problem is unlimited")
            return

        # Outgoing variable
        ratios = tableau[:-1, -1] / tableau[:-1, pivot_col]
        ratios[ratios <= 0] = np.inf  # Ignore negative 
        pivot_row = np.argmin(ratios)

        # Pivotage
        pivot_value = tableau[pivot_row, pivot_col]
        tableau[pivot_row] /= pivot_value  # Normalize pivot line

        for i in range(tableau.shape[0]):
            if i != pivot_row:
                tableau[i] -= tableau[i, pivot_col] * tableau[pivot_row]

        print_tableau(tableau)

    # Solution
    solution = np.zeros(num_vars)
    for i in range(num_constraints):
        if np.count_nonzero(tableau[i, :num_vars]) == 1:
            j = np.where(tableau[i, :num_vars] == 1)[0][0]
            solution[j] = tableau[i, -1]

    return solution, tableau[-1, -1]  # Solution and value of the objective function

def parse_args():
    parser = argparse.ArgumentParser(description="Solve a linear programming problem")
    
    # Arguments for the objective function
    parser.add_argument("--obj_x", type=float, required=True, help="Coefficient of x in the objective function")
    parser.add_argument("--obj_y", type=float, required=True, help="Coefficient of y in the objective function")
    parser.add_argument("--maximize", action="store_true", help="Maximize the function")
    
    # Arguments for constraints
    parser.add_argument("--constraint", action="append", nargs=3, metavar=('coeff_x', 'coeff_y','bound'),
                        help="Add a constraint like : coeff_x coeff_y bound", required=True)

    args = parser.parse_args()
    
    # Parsing of constraints
    constraints = []
    for constraint in args.constraint:
        coeff_x = float(constraint[0])
        coeff_y = float(constraint[1])
        bound = float(constraint[2])
        constraints.append((coeff_x, coeff_y, bound))
    
    # Coefficients of the objective function
    objective_coeffs = (args.obj_x, args.obj_y)
    
    return constraints, objective_coeffs, args.maximize

if __name__ == "__main__":
    constraints, objective_coeffs, maximize = parse_args()
    print('___USING OR TOOLS___\n')
    solve_with_ortools(constraints, objective_coeffs, maximize)
    print()

    print('___USING SIMPLEX STEP BY STEP___\n')
        
    num_constraints = len(constraints)
    num_vars = 2  # x et y
        
    A = np.zeros((num_constraints, num_vars))  
    b = np.zeros(num_constraints)  
    c = np.array([objective_coeffs[0], objective_coeffs[1]])  
        
    for i, (coeff_x, coeff_y, bound) in enumerate(constraints):
        A[i, 0] = coeff_x
        A[i, 1] = coeff_y
        b[i] = bound

    simplex_solution = simplex_step_by_step(c, A, b)
        
    if simplex_solution:
        print("\nSolution of the simplex methid (step by step) :")
        print("Solution : ", simplex_solution[0])
        print("Optimal value : ", simplex_solution[1])