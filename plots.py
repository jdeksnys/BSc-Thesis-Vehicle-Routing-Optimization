import matplotlib.pyplot as plt
import itertools
from TestCaseEnum import ROUTEPLOT
from matplotlib.lines import Line2D
from matplotlib.text import Text
import os




def PlotArcs(node_coords, main_arc_list, main_arc_dist, arc_list, plot_type, req_or_train_ids, orig_dict, dest_dict, duplicate_id, addStationLabels=True, addDistLabels=True, saveFigName=None, plot_vehicle_orig=False):
    """
    Constructs plot including nodes, node ids, distances, routes (vehicle or order).
    """
    plt.figure(figsize=(6, 6))

    for (m, n) in main_arc_list:
        edge = tuple(sorted([m, n]))
        linestyle = "-"
        x_values = [node_coords[m][0], node_coords[n][0]]
        y_values = [node_coords[m][1], node_coords[n][1]]
        x_mid_val = (node_coords[m][0] + node_coords[n][0]) / 2
        y_mid_val = (node_coords[m][1] + node_coords[n][1]) / 2
        plt.plot(x_values, y_values, color="black", linewidth=1, alpha=0.7, linestyle=linestyle)
        if(addDistLabels==True):
            plt.text(x_mid_val, y_mid_val, str(main_arc_dist[m,n]), color="grey")

    color_list = ["red", "blue", "green", "purple", "orange"]
    colors = itertools.cycle(color_list)

    edge_count = {}
    x_coords, y_coords = zip(*node_coords.values())

    if(addStationLabels==True):
        for key, (i, j) in node_coords.items():
            plt.text(i, j, str(key), fontsize=12, verticalalignment='bottom', horizontalalignment='right')

    def offset_line(x_values, y_values, offset=0.1):
        return [x + offset for x in x_values], [y + offset for y in y_values]

    plt.scatter(x_coords, y_coords, color='black', marker='o', label="Nodes")    

    for arcs, color in zip(arc_list, colors):
        for _, (m, n) in arcs.items():
            edge = tuple(sorted([m, n]))  # ensure (m, n) and (n, m) are treated the same
            
            if edge in edge_count:
                linestyle = "-"
                offset = 0.2 * edge_count[edge]
                x_values = [node_coords[m][0], node_coords[n][0]]
                y_values = [node_coords[m][1], node_coords[n][1]]
                x_values_offset, y_values_offset = offset_line(x_values, y_values, offset)
                plt.plot(x_values_offset, y_values_offset, color=color, linewidth=2, alpha=0.7, linestyle=linestyle)
            else:
                linestyle = "-"
                x_values = [node_coords[m][0], node_coords[n][0]]
                y_values = [node_coords[m][1], node_coords[n][1]]
                plt.plot(x_values, y_values, color=color, linewidth=2, alpha=0.7, linestyle=linestyle)
                edge_count[edge] = 1 
            
            edge_count[edge] += 1

    arc_legend = []
    arc_labels = []

    dummy_line = Line2D([0], [0], color='red', linewidth=0)  # line for red text (invisible)
    arc_legend.append(dummy_line)
    arc_labels.append("1: Node Id")
    arc_legend.append(dummy_line)
    arc_labels.append("2: Distance between nodes")

    line = Line2D([0], [0], color='black', linewidth=1)
    arc_legend.append(line)
    arc_labels.append(f"Existing arc")

    for i in range(0,len(req_or_train_ids)):
        c_id = i
        if(i>=len(color_list)):
            c_id = i-len(color_list)
        red_line = Line2D([0], [0], color=color_list[c_id], linewidth=2)
        arc_legend.append(red_line)

        if(plot_type == ROUTEPLOT.VEHICLE):
            vehicle = fr"$k={i}$"
            orig = fr"$o({i})$"
            
            if(orig_dict[req_or_train_ids[i]]>=duplicate_id):
                orig_val = orig_dict[req_or_train_ids[i]] - duplicate_id
            else:
                orig_val = orig_dict[req_or_train_ids[i]]

            route_label = fr"Vehicle {vehicle} route , {orig}={orig_val}"
        elif(plot_type == ROUTEPLOT.ORDER):
            order = fr"$r={i}$"
            orig = fr"$p({i})$"
            dest = fr"$d({i})$"
            route_label = fr"Order {order} route , ({orig}, {dest})=({orig_dict[req_or_train_ids[i]]}, {dest_dict[req_or_train_ids[i]]})"
        arc_labels.append(route_label)

    plt.xlabel("x coordinate")
    plt.ylabel("y coordinate")
    if(plot_type == ROUTEPLOT.VEHICLE):
        plt.title("Vehicle routes")
    elif(plot_type == ROUTEPLOT.ORDER):
         plt.title("Order routes")
    legend = plt.legend(arc_legend, arc_labels, loc='lower center', bbox_to_anchor=(0.5, 1.1), ncol=1)

    plt.subplots_adjust(top=0.6) 
    legend_bbox = legend.get_bbox_to_anchor()

    for text in legend.get_texts():
        if text.get_text() == "2: Distance between nodes":
            text.set_color('grey')

    plt.grid(True)


    if(saveFigName==None or saveFigName==""):
        plt.show()
    else:
        dir_name = os.path.dirname(saveFigName)
        base_name = os.path.splitext(os.path.basename(saveFigName))[0]
        ext = os.path.splitext(saveFigName)[1]
        i = 1
        new_base = base_name
        while any(os.path.splitext(f)[0] == new_base for f in os.listdir(dir_name)):
            new_base = f"{base_name} ({i})"
            i += 1

        new_path = os.path.join(dir_name, new_base + ext)
        plt.savefig(new_path)



from TestCaseEnum import TEST
import TestCaseAPI as API

def PlotGraphOnly(case, nodeIds=True, dist=False):
    N_coords = API.get_N_coords(case)
    N = API.get_N(N_coords)
    T = API.get_T(case)
    A, D = API.get_A_D(case)
    R, R_orig, R_dest, R_quant = [],[],[],[]
    K, K_orig, K_dest, K_cap, K_w_cost, K_l_cost = [],[],[],[],[],[]
    X_arc_list = []
    all_arcs_list = list(A)
    PlotArcs(N_coords, all_arcs_list, D, X_arc_list, ROUTEPLOT.VEHICLE, K, K_orig, K_dest, nodeIds, dist)