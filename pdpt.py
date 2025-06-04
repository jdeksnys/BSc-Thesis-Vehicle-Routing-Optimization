
import time
import gurobipy as gp
from gurobipy import GRB
import os
from TestCaseEnum import TEST, ROUTEPLOT
from PrintOrPlotEnum import PRINTORPLOT
import TestCaseAPI as API
from plots import PlotArcs, PlotGraphOnly
from parser import ParseStringsToDictsX, ParseStringsToDictsY
from itertools import chain, combinations
os.system('cls' if os.name == 'nt' else 'clear')




# EXECUTION PARAMETERS /////////////////////////////////////////////////////////////
"""
Set output type: PRINT, PLOT, PLOTSAVE, PLOTGRAPHONLY.
Select test case (differ by data: nodes, coordinates, vehicles, orders).
"""
output_type = PRINTORPLOT.PLOT
case = TEST.CASE_8




# GUROBI /////////////////////////////////////////////////////////////
"""
set Gurobi solver parameters to use basic simplex method
"""
start_time = time.time()
m = gp.Model("pdp-transshipments")
m.setParam(GRB.Param.LazyConstraints, 1) # for subtour elimination using subsets of used nodes

m.setParam(GRB.Param.Presolve, 0)       # off
m.setParam(GRB.Param.Heuristics, 0)     # off
m.setParam(GRB.Param.Cuts, 0)           # off
m.setParam(GRB.Param.Method, 0)         # primal simplex (regular simplex described in )




# PARAMETERS /////////////////////////////////////////////////////////
"""
Get parameter values, which depend on the test case.
Set omega value which is used for identifying artificial nodes for the workaroud (allowing to pass origin node more than once).
"""

if(output_type == PRINTORPLOT.PLOTGRAPHONLY):
    PlotGraphOnly(case)
    exit()

N_coords = API.get_N_coords(case)
N = API.get_N(N_coords)
T = API.get_T(case)
A, D = API.get_A_D(case)
R, R_orig, R_dest, R_quant, R_w_cost = API.get_R_R_orig_R_dest_R_quant_R_w_cost(case)
K, K_orig, K_dest, K_cap, K_l_cost = API.get_K_K_orig_K_dest_K_cap_K_l_cost(case, False)
duplicate_id = 1000 # omega value for allowing to pass origin node more than once
K_orig, N_coords, N, A, D = API.DoubleTrainOriginNodes(K_orig, N_coords, N, A, D, duplicate_id)




# DECISION VARIABLES /////////////////////////////////////////////////
"""
Add decision variables and their data type.
"""
# 1 if train k uses arc a. 0 otherwise.
X = m.addVars(K, A, vtype=GRB.BINARY, name="X")


# 1 if request r transported by tran k uses arc a. 0 otherwise.
Y = m.addVars(K, R, A, vtype=GRB.BINARY, name="Y")

# 1 if request r is transferred from train k to train l!=k at node i. 0 otherwise.
S = m.addVars(K, K, N, R, vtype=GRB.BINARY, name="S")



# OBJECTIVE FUNCTION /////////////////////////////////////////////////
"""
Set objective function, which consists of two summations (total wagon travel cost and total locomotive travel cost).
"""
objective_expr1 = gp.LinExpr()
objective_expr2 = gp.LinExpr()

for k in K:
    for a in A:
        for r in R:
            objective_expr1 += R_w_cost[r] * D[a] * Y[k, r, a[0], a[1]] * R_quant[r]
for k in K:
    for a in A:
        objective_expr2 += K_l_cost[k] * D[a] * X[k,a[0],a[1]]

objective_expr = objective_expr1 + objective_expr2
m.setObjective(objective_expr, GRB.MINIMIZE)



# CONSTRAINTS ////////////////////////////////////////////////////////
"""
Define mathematical model constraints, except constraints for ensuring decision variables are binary.
"""

# constr_5
for k in K:
    i = K_orig[k]
    m.addConstr(
        gp.quicksum(
            X.select(k,i,"*")
        ) <= 1
    )

# constr_6
for r in R:
    i = R_orig[r]
    m.addConstr(
        gp.quicksum(
            gp.quicksum(
                Y.select(k,r,i,"*")
            ) for k in K
        ) == 1
    )

# constr_7
for r in R:
    i = R_dest[r]
    m.addConstr(
        gp.quicksum(
            gp.quicksum(
                Y.select(k,r,"*",i)
            ) for k in K
        ) == 1
    )

# constr_8
for k in K:
    N_filtered = [n for n in N if n!=K_orig[k] and n!=K_dest[k]]
    for i in N_filtered:
        m.addConstr(
            gp.quicksum(
                X.select(k, i, "*")
            ) - gp.quicksum(
                X.select(k, "*", i)
            ) == 0
        )

# constr_9
for r in R:
    T_filtered = [t for t in T if t!=R_orig[r] and t!= R_dest[r]]
    for i in T_filtered:
        m.addConstr(
            gp.quicksum(
                gp.quicksum(
                    Y.select(k,r,i,"*")
                ) for k in K
            ) - gp.quicksum(
                gp.quicksum(
                    Y.select(k,r,"*",i)
                ) for k in K
            ) == 0
        )

# constr_10
for k in K:
    for r in R:
        N_filtered = [n for n in N if n not in T and n!=R_orig[r] and n!=R_dest[r]]
        for i in N_filtered:
            m.addConstr(
                gp.quicksum(
                    Y.select(k,r,i,"*")
                ) - gp.quicksum(
                    Y.select(k,r,"*",i)
                ) == 0
            )

# constr_11
for a in A:
    for k in K:
        for r in R:
            m.addConstr(
                Y[k,r,a[0],a[1]] <= X[k,a[0],a[1]]
            )

# constr_12
for r in R:
    for i in T:
        for k in K:
            for l in K:
                if(l==k):
                    continue
                m.addConstr(
                    gp.quicksum(
                        Y.select(k,r,"*",i)
                    ) + gp.quicksum(
                        Y.select(l,r,i,"*")
                    ) <= S[k,l,i,r] + 1
                )

# constr_13
for r in R:
    for i in T:
        for k in K:
            for l in K:
                if(k == l):
                    continue
                m.addConstr(
                    S[k,l,i,r] <= gp.quicksum(
                        Y.select(k,r,"*",i)
                    )
                )

# constr_14
for r in R:
    for i in T:
        for k in K:
            for l in K:
                if(k == l):
                    continue
                m.addConstr(
                    S[k,l,i,r] <= gp.quicksum(
                        Y.select(l,r,i,"*")
                    )
                )

# constr_15
for a in A:
    for k in K:
        m.addConstr(
            gp.quicksum(
                R_quant[r] * Y.select(k,"*",a[0],a[1])
            )
            <= K_cap[k] * X[k,a[0],a[1]]
        )


# constr_16
for k in K:
    i = K_orig[k]
    m.addConstr(
        gp.quicksum(
            X.select(k,"*","*")
        ) <= (
            # (indicator if vehicle goes: duplicate_origin --> origin) * ()
            X[k,i,i-duplicate_id] * len(N)
        )
    )

 


# constr_17, constr_18
"""
Subtour elimination callbacks.
"""
class SubtourElimCallback:
    def __init__(self, X, Y, S):
        self.X = X
        self.Y = Y
        self.S = S
    
    def __call__(self, model, where):
        if where == GRB.Callback.MIPSOL:
            try:
                self.eliminate_subtours(model)
            except Exception:
                print(f"\n\n\n[ERROR] {Exception}")
                model.terminate()
        
    def eliminate_subtours(self, model):
        x_sol = {key: model.cbGetSolution(var) for key, var in self.X.items()}
        y_sol = {key: model.cbGetSolution(var) for key, var in self.Y.items()}
        used_trains = set([key[0] for key,val in x_sol.items() if val==1])
        used_trains = used_trains.union(set([key[0] for key,val in y_sol.items() if val==1]))

        def all_subsets(s):
            return list(chain.from_iterable(combinations(s, r) for r in range(1, len(s))))

        for k in used_trains:
            used_nodes_x = set(
                [key[1] for key,val in x_sol.items() if val==1 and key[0]==k] 
                + 
                [key[2] for key,val in x_sol.items() if val==1 and key[0]==k]
            )
            used_node_subsets_x = all_subsets(used_nodes_x)
            for subset in used_node_subsets_x:
                arcs_outside_subset = (gp.quicksum(
                        x_sol[k, i, j]
                        for (i,j) in A
                        if i in subset and j not in subset
                    ) + gp.quicksum(
                        x_sol[k, i, j]
                        for (i,j) in A
                        if j in subset and i not in subset
                    )).getValue()
                if(arcs_outside_subset < 0.5):
                    expr = gp.quicksum(
                        self.X[k, i, j]
                        for (i,j) in A
                        if i in subset and j not in subset
                    ) + gp.quicksum(
                        self.X[k, i, j]
                        for (i,j) in A
                        if j in subset and i not in subset
                    )
                    model.cbLazy(expr >= 1)

            requests_for_k = set([key[1] for key,val in y_sol.items() if val==1 and key[0]==k])
            for r in requests_for_k:
                used_nodes_y_for_k_for_r = set(
                    [key[2] for key,val in y_sol.items() if val==1 and key[0]==k and key[1]==r] 
                    + 
                    [key[3] for key,val in y_sol.items() if val==1 and key[0]==k and key[1]==r]
                )
                used_node_subsets_y = all_subsets(used_nodes_y_for_k_for_r)
                for subset in used_node_subsets_y:
                    arcs_outside_subset = (gp.quicksum(
                            y_sol[k,r,i,j]
                            for (i,j) in A
                            if i in subset and j not in subset
                        ) + gp.quicksum(
                            y_sol[k,r,i,j]
                            for (i,j) in A
                            if j in subset and i not in subset
                        )).getValue()
                    if(arcs_outside_subset < 0.5):
                        expr = gp.quicksum(
                            self.Y[k,r,i,j]
                            for (i,j) in A
                            if i in subset and j not in subset
                        ) + gp.quicksum(
                            self.Y[k,r,i,j]
                            for (i,j) in A
                            if j in subset and i not in subset
                        )
                        model.cbLazy(expr >= 1)





# RUN ////////////////////////////////////////////////////////////////
"""
Specify callback functions (some constraints are implemented as callbacks).
Execute the Gurobi solver.
"""
callback = SubtourElimCallback(X, Y, S)
m.update()

auto_gen_plot_dir = os.path.join(os.getcwd(), 'auto_gen_plots', 'model.lp')
os.makedirs(auto_gen_plot_dir, exist_ok=True)
auto_gen_plot_path = os.path.join(auto_gen_plot_dir, 'model.lp')

m.write(auto_gen_plot_path)  # or "model.mps"
m.optimize(callback)

end_time = time.time()
print(f"---------------Runtime: {end_time - start_time:.6f} seconds---------------")





# PRINT OR PLOT COMMANDS ////////////////////////////////////////////////////////////////
"""
Print or plot results.
"""
if(m.Status == GRB.OPTIMAL):
    num_solutions = m.SolCount
    vars = m.getVars()
    for i in range(num_solutions):
        m.setParam(GRB.Param.SolutionNumber, i)
        decisionVars = m.getVars()
        X_arc_list_str = []
        Y_arc_list_str = []
        for v in decisionVars:
            if(v.x != 0):
                if(output_type == PRINTORPLOT.PRINT):
                    if(v.VarName[0]=="S"):
                        continue
                    elif(v.VarName[0]=="Z"):
                        continue
                    print(f'{v.varName}: {int(v.x)}')
                elif(output_type == PRINTORPLOT.PLOT or output_type == PRINTORPLOT.PLOTSAVE):
                    if(v.VarName[0]=="X"):
                        X_arc_list_str.append(v.varName)
                    elif(v.VarName[0]=="Y"):
                        Y_arc_list_str.append(v.varName)
        if(output_type == PRINTORPLOT.PLOT or output_type == PRINTORPLOT.PLOTSAVE):
            if(output_type == PRINTORPLOT.PLOT):
                fig1_path = None
                fig2_path = None
            else:
                fig1_path= os.path.join(os.path.expanduser("~"), "Desktop") + "/auto_gen_plots/output_type_t.png"
                fig2_path= os.path.join(os.path.expanduser("~"), "Desktop") + "/auto_gen_plots/output_type_w.png"
            X_arc_list = ParseStringsToDictsX(X_arc_list_str)
            Y_arc_list = ParseStringsToDictsY(Y_arc_list_str)
            all_arcs_list = list(A)
            
            PlotArcs(N_coords, all_arcs_list, D, X_arc_list, ROUTEPLOT.VEHICLE, K, K_orig, K_dest, duplicate_id, True, True, saveFigName=fig1_path, plot_vehicle_orig=True)
            PlotArcs(N_coords, all_arcs_list, D, Y_arc_list, ROUTEPLOT.ORDER, R, R_orig, R_dest, duplicate_id, True, True, saveFigName=fig2_path, plot_vehicle_orig=True)
else:
    print("FAIL")

