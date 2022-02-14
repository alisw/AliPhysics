import uproot
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm


shift_list = [[1, -2], [2, -1], [3, 1], [4, 2]]
dcaxy_list = [[0, 1.0], [1, 1.4]]
dcaz_list = [[0, 0.5], [1, 0.75], [2, 1.25], [3, 1.50]]
pid_list = [[0, 3.25], [1, 3.5]]
tpc_list = [[0, 60], [1, 65], [2, 75], [3, 80]]
width_list = [[4, -2], [2, -1], [1, +1], [3, 2]]
cuts = {"shift": [shift_list, "Bin shift"], "dcaxy": [dcaxy_list, "$DCA_{xy}$ (mm)"], "dcaz": [dcaz_list, "$DCA_{z}$ (cm)"], "pid": [
    pid_list, "$n\sigma_{TPC}$"], "tpc": [tpc_list, "TPC clusters"], "width": [width_list, "Bin width"]}

inFile = uproot.open("spectra.root")
normHist = inFile["nuclei_deuterons_/deuterons/9/Joined/JoinedSpectraM9"]

norm = normHist.values[:-3]
pt = [0.5 * (x + y) for x, y in zip(normHist.edges[:-4], normHist.edges[1:-3])]


colors = cm.rainbow(np.linspace(0, 1, len(norm)))
cMap = plt.get_cmap('jet')
cNorm = matplotlib.colors.Normalize(vmin=min(pt), vmax=max(pt))
scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cMap)

for key, record in cuts.items():
    x = np.array([])
    y = np.array([])
    fig,ax = plt.subplots()
    for obj_list in record[0]:
        obj_label=obj_list[0]
        obj_val=obj_list[1]
        values = inFile["nuclei_deuterons_{}{}/deuterons/9/Joined/JoinedSpectraM9".format(
            key, obj_label)].values[:-3]
        values = values / norm
        x = np.append(x, np.array([obj_val for _ in range(0, len(values))]))
        y = np.append(y, values)
        plt.scatter(np.array([obj_val for _ in range(0, len(values))]),
                    values, color=scalarMap.to_rgba(pt), edgecolors='none')

    scalarMap.set_array(pt)
    fig.colorbar(scalarMap).set_label("$p_{T}$ (GeV/$c$)")
    plt.ylabel("$S_{var}$ / $S_{ref}$")
    plt.xlabel(record[1])
    ax.tick_params(axis="y",direction="in")
    ax.tick_params(axis="x",direction="in")
    ax.xaxis.set_label_coords(0.9,-0.07)
    ax.yaxis.set_label_coords(-0.115,0.9)
    plt.text(0.5, 0.92, 'This work', ha='center', va='center', transform=ax.transAxes, fontweight='bold', fontsize=14) 
    
    print(np.corrcoef(x, y))

    plt.savefig("{}.pdf".format(key))


# for name, keys in cuts:
