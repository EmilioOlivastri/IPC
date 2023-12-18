import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

from texttable import Texttable
import latextable

def latexfy(table_array, dataset):
    # Example 3 - Position
    table_3 = Texttable()
    table_3.set_cols_align(["c"] * 11)
    table_3.set_deco(Texttable.HEADER | Texttable.VLINES)
    table_3.add_rows(table_array)
    print('Texttable Output:')
    print(table_3.draw())
    print('\nLatextable Output:')
    print(latextable.draw_latex(table_3, caption=dataset, label="table:position", position='ht'))

def read_PRT(path_pr) :

    data = None

    with open(path_pr) as f:
        lines = f.readlines()
        pr, rec = lines[0].split()
        data = [float(pr), float(rec), float(lines[1][:-1])]

    return np.array(data)

def read_tj_info(path_tf_info) :

    data = None

    with open(path_tf_info) as f:
        lines = f.readlines()
        _, ate, _, abs_t, _, rpe, _ = lines[0].split()
        data = [float(ate), float(abs_t), float(rpe)]

    return np.array(data)


def getData(base_path, datasets, outliers, runs, alg_gtsam, alg_g2o, alg_rpgo) :

    final_dict = {'ALG' : [],
                  'DSET': [],
                  'OUTLIERS' :[],
                  'RUN' :[],
                  'PREC' :[],
                  'REC' :[],
                  'F1': [],
                  'TIME' :[],
                  'ATE' : [],
                  'ABS_T': [],
                  'RPE' : []}
    exten = ['.TE', '.PR']

    for dataset in datasets :
        first_sec = base_path + dataset + '/RESULTS/'
        for alg in alg_gtsam :
            second_sec = first_sec + 'GTSAM_' + alg + '/'
            for out in outliers :
                third_sec = second_sec + out + '/'
                for run in runs :
                    path_prt = third_sec + run + exten[1]
                    path_tf_info = third_sec + run + exten[0]
                    res_prt = read_PRT(path_prt)
                    res_tf = read_tj_info(path_tf_info)
                    final_dict['ALG'].append(alg) 
                    final_dict['DSET'].append(dataset) 
                    final_dict['OUTLIERS'].append(int(out))    
                    final_dict['RUN'].append(int(run))    
                    final_dict['PREC'].append(float(res_prt[0]))    
                    final_dict['REC'].append(float(res_prt[1]))
                    f1 = (2 * float(res_prt[0]) * float(res_prt[1])) / (float(res_prt[0]) + float(res_prt[1]))
                    final_dict['F1'].append(f1)   
                    final_dict['TIME'].append(float(res_prt[2]))    
                    final_dict['ATE'].append(float(res_tf[0]))
                    final_dict['ABS_T'].append(float(res_tf[1]))
                    final_dict['RPE'].append(float(res_tf[2])) 

        for alg in alg_g2o :
            second_sec = first_sec + 'G2O_' + alg + '/'
            for out in outliers :
                third_sec = second_sec + out + '/'
                for run in runs :
                    path_prt = third_sec + run + exten[1]
                    path_tf_info = third_sec + run + exten[0]
                    res_prt = read_PRT(path_prt)
                    res_tf = read_tj_info(path_tf_info)
                    final_dict['ALG'].append(alg)
                    final_dict['DSET'].append(dataset) 
                    final_dict['OUTLIERS'].append(int(out))    
                    final_dict['RUN'].append(int(run))    
                    final_dict['PREC'].append(float(res_prt[0]))    
                    final_dict['REC'].append(float(res_prt[1]))
                    f1 = (2 * float(res_prt[0]) * float(res_prt[1])) / (float(res_prt[0]) + float(res_prt[1]))
                    final_dict['F1'].append(f1)
                    final_dict['TIME'].append(float(res_prt[2]))    
                    final_dict['ATE'].append(float(res_tf[0]))
                    final_dict['ABS_T'].append(float(res_tf[1]))
                    final_dict['RPE'].append(float(res_tf[2])) 

        for alg in alg_rpgo :
            second_sec = first_sec + 'RPGO_' + alg + '/'
            for out in outliers :
                third_sec = second_sec + out + '/'
                for run in runs :
                    path_prt = third_sec + run + exten[1]
                    path_tf_info = third_sec + run + exten[0]
                    res_prt = read_PRT(path_prt)
                    res_tf = read_tj_info(path_tf_info)
                    final_dict['ALG'].append(alg)
                    final_dict['DSET'].append(dataset) 
                    final_dict['OUTLIERS'].append(int(out))    
                    final_dict['RUN'].append(int(run))    
                    final_dict['PREC'].append(float(res_prt[0]))    
                    final_dict['REC'].append(float(res_prt[1]))
                    f1 = (2 * float(res_prt[0]) * float(res_prt[1])) / (float(res_prt[0]) + float(res_prt[1]))
                    final_dict['F1'].append(f1)
                    final_dict['TIME'].append(float(res_prt[2]))    
                    final_dict['ATE'].append(float(res_tf[0]))
                    final_dict['ABS_T'].append(float(res_tf[1]))
                    final_dict['RPE'].append(float(res_tf[2])) 


    return pd.DataFrame.from_dict(final_dict)

def plot(results, vot_results, xlabel, col, methods_name, colors, markers, terror_min = -1.0, terror_max = 100,	rerror_min = 0.01,	rerror_max = 45):
    
    plt.subplot(2, 2, 1); 
    for idx, alg in enumerate(methods_name):
        ticks = []
        means = []
        for data in results[alg] :
            ticks.append(data['out'])
            means.append(data['ate'][0])
        
        plt.plot(ticks, means, label=alg, marker=markers[idx], color=colors[idx], markersize=8)
    
    #plt.ylim(ymin=terror_min, ymax=terror_max)
    #plt.xlabel(xlabel)
    plt.grid()
    
    plt.subplot(2, 2, 2)
    for idx, alg in enumerate(methods_name):
        ticks = []
        means = []
        for data in results[alg] :
            ticks.append(data['out'])
            means.append(data['rpe'][0])
        
        plt.plot(ticks, means, label=alg, marker=markers[idx], color=colors[idx], markersize=8)
    
    #plt.ylim(ymin=rerror_min, ymax=rerror_max)	
    #plt.xlabel(xlabel)
    plt.grid()

    plt.subplot(2, 1, 2)
    idx_alg = 0
    for idx, alg in enumerate(methods_name):
        ticks = []
        means = []
        if alg == 'VOT' :
            idx_alg = idx
            continue
        for data in results[alg] :
            ticks.append(data['out'])
            means.append(data['time'][0])
        
        plt.plot(ticks, means, label=alg, marker=markers[idx], color=colors[idx], markersize=8)

    ticks_vt = [str(out * (10) + 10) for out in range(10)]
    means_vt = []
    for out in range(10):
        tmp = []
        for dataset in vot_results :
            tmp.append(vot_results[dataset][out]['total_time'][0])
        means_vt.append(np.mean(tmp))

    plt.plot(ticks_vt, means_vt, label="VOT", marker=markers[idx_alg], color=colors[idx_alg], markersize=8)
    #print(ticks_vt)
    #print(means_vt)

    ##plt.ylim(ymin=rerror_min, ymax=rerror_max)
    plt.xlabel(xlabel)	
    plt.grid()


def printMethodStats(diz, method):

    elems = [1, 4, 6, 9]
    for x in elems :
        print(f'out = {diz[method][x]["out"][0]} | f1 = {diz[method][x]["f1"][0]} | prec = {diz[method][x]["prec"][0]} | rpe = {diz[method][x]["rpe"][0]}')
        print("-----------")
    
def main() :

    base_path = '/home/slam-emix/Datasets/ICRA_RESULTS/'
    outliers = ['10', '20', '30', '40', '50', '60', '70', '80', '90', '100']
    runs = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']
    datasets = ['MIT', 'INTEL', 'M3500', 'CSAIL', 'FRH', 'FR079']
    dataset_inliers = [20.0, 256.0, 1954.0, 128.0, 229.0, 1505.0]
    opt = ['G2O', 'GTSAM']
    alg_gtsam = ['GNC', 'GM', 'HUBER', 'DCS']
    #alg_gtsam = ['GNC', 'HUBER', 'DCS']
    #alg_g2o = ['VOT', 'MAXMIX']
    alg_g2o = ['MAXMIX', 'VOT', 'ADAPT']
    alg_rpgo = ['PCM']
    full_list_algs = alg_gtsam + alg_g2o + alg_rpgo
    markers = ['X', 'o', '^', 'D', 'v', 's', 'p', 'P']
    colors = ['red', 'blue', 'green', 'purple', 'pink', 'cyan', 'orange', 'brown']
                
    df = getData(base_path, datasets, outliers, runs, alg_gtsam, alg_g2o, alg_rpgo)
    #print(df.head())

    diz = {}
    for alg in full_list_algs :
        alg_df = df[df["ALG"] == alg]
        alg_df = alg_df[alg_df["DSET"] == "FR079"]
        diz[alg] = []

        # compute statistics for elements of interest
        for out in outliers :
            expr_df = alg_df[alg_df["OUTLIERS"] == int(out)]
            diz[alg].append({'out': out, 
                            'prec' : [expr_df["PREC"].mean(), expr_df["PREC"].std()],
                            'rec' : [expr_df["REC"].mean(), expr_df["REC"].std()],
                            'f1' : [expr_df["F1"].mean(), expr_df["F1"].std()],
                            'ate' : [expr_df["ATE"].mean(), expr_df["ATE"].std()],
                            'rpe' : [expr_df["RPE"].mean(), expr_df["RPE"].std()],
                            'time': [expr_df["TIME"].mean(), expr_df["TIME"].std()]})


    printMethodStats(diz, "MAXMIX")
    print("#########")
    printMethodStats(diz, "VOT")

    vot_df = df[df["ALG"] == "VOT"] 
    vot_diz = {}
    for idd, dataset in enumerate(datasets):
        dts_df = vot_df[vot_df["DSET"] == dataset]
        vot_diz[dataset] = []
        increment = 1.1
    
        # compute statistics for elements of interest
        for out in outliers :
            total_edges = int(increment * dataset_inliers[idd])
            expr_df = dts_df[dts_df["OUTLIERS"] == int(out)]
            mean_vot = expr_df["TIME"].mean()
            std_vot = expr_df["TIME"].std()
            vot_diz[dataset].append({'out': out, 
                                     'prec' : [expr_df["PREC"].mean(), expr_df["PREC"].std()],
                                     'rec' : [expr_df["REC"].mean(), expr_df["REC"].std()],
                                     'f1' : [expr_df["F1"].mean(), expr_df["F1"].std()],
                                     'rpe' : [expr_df["RPE"].mean(), expr_df["RPE"].std()],
                                     'voting_time' : [mean_vot, std_vot],
                                     'total_time': [total_edges * mean_vot, std_vot]})
            increment += 0.1
    #print(vot_diz['MIT'])     

    # Computing table
    table_dt, table_tt  = [], []
    table_dt.append(['Dataset',  '10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%', '90\%', '100\%'])
    table_tt.append(['Dataset',  '10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%', '90\%', '100\%'])
    for idd, dataset in enumerate(datasets):
            row_dt, row_tt = [dataset], [dataset]
            for data in vot_diz[dataset] :
                row_dt.append(data['voting_time'][0])
                row_tt.append(data['total_time'][0])

            table_dt.append(row_dt)
            table_tt.append(row_tt)

    #latexfy(table_dt,'Time')
    #latexfy(table_tt,'F1')
    
    fig = plt.figure(figsize=(12, 4))
    plt.subplot(2, 2, 1)
    plt.ylabel('ATE [m]')
    plt.subplot(2, 2, 2)
    plt.ylabel('RPE [m]')
    plt.subplot(2, 1, 2)
    plt.ylabel('Convergence [s]')


    plot(results=diz, vot_results=vot_diz, xlabel='Outliers [%]', col=1, methods_name=full_list_algs, 
         colors=colors, markers=markers, terror_min= 0.45, terror_max=1.0, 
         rerror_min=0.3, rerror_max=1.0)

    plt.tight_layout()
    #plt.show()

    return 



if __name__ == '__main__' :
    main()