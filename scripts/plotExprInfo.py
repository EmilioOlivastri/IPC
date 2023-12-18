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


def classic_boxplot(dataframe, key_xyh, ax, algs, colors, y_coord) :

    # Computing symbol positioning for boxplot
    w = 0.5
    central_ticks = np.arange(0.0, 10, step=1.0)
    l_list = len(algs)
    y = np.ones(len(central_ticks)) * y_coord
    central_elem = []
    left_elem = []
    right_elem = []
    palette_c = []
    palette_l = []
    palette_r = []
    if ( l_list % 2 == 1 ) :
        c_idx = int(l_list / 2)
        central_elem.append(algs[c_idx])
        palette_c.append(colors[c_idx])

      
    for idx in range(int(l_list/2)) :
        right_elem.append(algs[idx])
        palette_r.append(colors[idx])
        left_elem = [algs[-idx -1]] + left_elem
        palette_l = [colors[-idx - 1]] + palette_l

    real_label_order = left_elem + central_elem + right_elem
    palette = palette_l + palette_c + palette_r

    # Symbol positioning completed

    box_plt = sns.boxplot(x=dataframe[key_xyh[0]], y=dataframe[key_xyh[1]], hue=dataframe[key_xyh[2]], 
                          hue_order=real_label_order, palette=palette, width=1.5, fliersize=0)
    #plt.yscale('symlog')
    ax.set_xticks(np.arange(-0.5, 10, step=1.0))
    ax.set_xticks(central_ticks, minor=True, labels=['10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    plt.grid()
    ax.legend() 

    return


def simple_plot(data_diz, algs, colors, markers,name, type_d) :

    fig = plt.figure(name)
    idx = 0
    for alg in algs :
        ticks = []
        means = []
        stds = []
        for data in data_diz[alg] :
            ticks.append(data['out'])
            means.append(data[type_d][0])
            stds.append(data[type_d][1])

        plt.plot(ticks, means, linewidth=2.5, color=colors[idx], label=alg, marker=markers[idx],markersize=12)
        idx += 1

    plt.xlabel('%outliers', fontsize=16)
    plt.ylabel(type_d, fontsize=16)
    #plt.yticks([0.25, 0.5, 0.75, 1.0], [0.25, 0.5, 0.75, 1.0])
    plt.grid()
    return


def daltonic_boxplot(dataframe, key_xyh, ax, algs, colors, markers, y_coord, turn_off=False) :

    # Computing symbol positioning for boxplot
    w = 0.5
    central_ticks = np.arange(0.0, 10, step=1.0)
    l_list = len(algs)
    step = w / l_list
    starting_off = 0.5 * step
    y = np.ones(len(central_ticks)) * y_coord
    central_elem = []
    left_elem = []
    right_elem = []
    palette_c = []
    palette_l = []
    palette_r = []
    if ( l_list % 2 == 1 ) :
        starting_off = step
        c_idx = int(l_list / 2)
        ax.scatter(central_ticks, y, c=colors[c_idx], 
                        marker=markers[c_idx], s=60, label=algs[c_idx])
        central_elem.append(algs[c_idx])
        palette_c.append(colors[c_idx])

      
    for idx in range(int(l_list/2)) :
        xpos1 = central_ticks + starting_off + idx * step
        xpos2 = central_ticks - starting_off - idx * step 
        ax.scatter(xpos1, y, c=colors[idx], 
                        marker=markers[idx], s=60, label=algs[idx])
        ax.scatter(xpos2, y, c=colors[-idx - 1], 
                        marker=markers[-idx- 1], s=60, label=algs[-idx -1])
        right_elem.append(algs[idx])
        palette_r.append(colors[idx])
        left_elem = [algs[-idx -1]] + left_elem
        palette_l = [colors[-idx - 1]] + palette_l

    real_label_order = left_elem + central_elem + right_elem
    palette = palette_l + palette_c + palette_r

    # Symbol positioning completed

    box_plt = sns.boxplot(x=dataframe[key_xyh[0]], y=dataframe[key_xyh[1]], hue=dataframe[key_xyh[2]], 
                          hue_order=real_label_order, palette=palette, width=w)
    #plt.yscale('symlog')
    ax.set_xticks(np.arange(-0.5, 10, step=1.0))
    ax.set_xticks(central_ticks, minor=True, labels=['10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    if turn_off == False :
        plt.yticks(np.arange(0.0, 1.25, step=0.25))
    plt.grid()
    h, l = ax.get_legend_handles_labels()
    ax.legend(h[:l_list], l[:l_list]) 

    return 

def getData(base_path, datasets, outliers, runs, alg_gtsam, alg_g2o, alg_rpgo) :

    final_dict = {'ALG' : [], 
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

def plot(results, xlabel, col, methods_name, colors, markers, terror_min = -1.0, terror_max = 100,	rerror_min = 0.01,	rerror_max = 45):
    
    plt.subplot(2, 2, 1); 
    for idx, alg in enumerate(methods_name):
        ticks = []
        means = []
        for data in results[alg] :
            ticks.append(data['out'])
            means.append(data['prec'][0])
        
        plt.plot(ticks, means, label=alg, marker=markers[idx], color=colors[idx], markersize=8)
    
    plt.ylim(ymin=terror_min, ymax=terror_max)
    plt.grid()
    
    plt.subplot(2, 2, 2)
    for idx, alg in enumerate(methods_name):
        ticks = []
        means = []
        for data in results[alg] :
            ticks.append(data['out'])
            means.append(data['rec'][0])
        
        plt.plot(ticks, means, label=alg, marker=markers[idx], color=colors[idx], markersize=8)
    
    plt.ylim(ymin=rerror_min, ymax=rerror_max)	
    plt.grid()
    #plt.xlabel(xlabel)

    plt.subplot(2, 1, 2)
    for idx, alg in enumerate(methods_name):
        ticks = []
        means = []
        for data in results[alg] :
            ticks.append(data['out'])
            means.append(data['f1'][0])
        
        plt.plot(ticks, means, label=alg, marker=markers[idx], color=colors[idx], markersize=8)

    plt.xlabel(xlabel)
    ##plt.ylim(ymin=rerror_min, ymax=rerror_max)	
    plt.grid()


def main() :

    base_path = '/home/slam-emix/Datasets/ICRA_RESULTS/'
    outliers = ['10', '20', '30', '40', '50', '60', '70', '80', '90', '100']
    runs = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']
    datasets = ['MIT', 'INTEL', 'M3500', 'CSAIL', 'FRH', 'FR079']
    opt = ['G2O', 'GTSAM']
    #alg_gtsam = ['GNC', 'GM', 'HUBER', 'DCS']
    alg_gtsam = ['GNC', 'HUBER', 'DCS']
    alg_g2o = ['VOT', 'MAXMIX']
    #alg_g2o = ['MAXMIX', 'VOT', 'ADAPT']
    alg_rpgo = ['PCM']
    full_list_algs = alg_gtsam + alg_g2o + alg_rpgo
    markers = ['X', 'o', '^', 'D', 'v', 's', 'p', 'P']
    colors = ['red', 'blue', 'green', 'purple', 'pink', 'cyan', 'orange', 'brown']
                
    df = getData(base_path, datasets, outliers, runs, alg_gtsam, alg_g2o, alg_rpgo)
    #print(df.head())

    diz = {}
    for alg in full_list_algs :
        alg_df = df[df["ALG"] == alg]
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

    #print(diz['GNC'])     

    # Computing table
    table_p, table_r, table_dt, table_f1  = [], [], [], []
    table_p.append(['Algorithm',  '10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%', '90\%', '100\%'])
    table_r.append(['Algorithm',  '10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%', '90\%', '100\%'])
    table_dt.append(['Algorithm', '10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%', '90\%', '100\%'])
    table_f1.append(['Algorithm', '10\%', '20\%', '30\%', '40\%', '50\%', '60\%', '70\%', '80\%', '90\%', '100\%'])
    for alg in full_list_algs:
            row_p, row_r, row_t, row_f1 = [alg], [alg], [alg], [alg]
            for data in diz[alg] :
                row_p.append(data['prec'][0])
                row_r.append(data['rec'][0])
                row_t.append(data['time'][0])
                row_f1.append(data['rpe'][0])

            table_p.append(row_p)
            table_r.append(row_r)
            table_dt.append(row_t)
            table_f1.append(row_f1)

    latexfy(table_p, 'Prec')
    latexfy(table_r, 'Rec')
    latexfy(table_dt,'Time')
    latexfy(table_f1,'F1')

    #simple_plot(diz, full_list_algs, colors, markers,'ATE [m]', 'ate')
    #simple_plot(diz, full_list_algs, colors, markers,'RPE [M]', 'rpe')
    
    fig = plt.figure(figsize=(12, 4))
    plt.subplot(2, 2, 1)
    plt.ylabel('Precision')
    plt.subplot(2, 2, 2)
    plt.ylabel('Recall')
    plt.subplot(2, 1, 2)
    plt.ylabel('F1')

    #plt.subplot(3, 1, 1)
    #plt.ylabel('Time [s]')
    #plt.subplot(3, 1, 2)
    #plt.ylabel('ATE [m]')
    #plt.subplot(3, 1, 3)
    #plt.ylabel('RPE [m]')

    plot(results=diz, xlabel='Outliers [%]', col=1, methods_name=full_list_algs, 
         colors=colors, markers=markers, terror_min= 0.45, terror_max=1.0, 
         rerror_min=0.3, rerror_max=1.0)

    ax = plt.subplot(2, 2, 1)
    handles, labels = ax.get_legend_handles_labels()
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.25)
    fig.legend(handles, labels, loc='lower center', ncol=3)

    plt.show()


    #print(diz['GNC'].head())
    #print(diz['GNC'].mean())
    #print(diz['GNC'].std())


    ## BOX PLOTS ##

    # PREC #
    #fig = plt.figure('PrecisionBox')
    #ax_prec = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    ##daltonic_boxplot(df, ['OUTLIERS', 'PREC', 'ALG'], ax_prec, full_list_algs, colors, markers, 1.01)
    #classic_boxplot(df, ['OUTLIERS', 'PREC', 'ALG'], ax_prec, full_list_algs, colors, 1.0)
    
    '''
    fig_rec = plt.figure('Recall')
    ax_rec = fig_rec.add_axes([0.1, 0.1, 0.8, 0.8])
    #daltonic_boxplot(df, ['OUTLIERS', 'REC', 'ALG'], ax_rec, full_list_algs, colors, markers, 1.005)
    classic_boxplot(df, ['OUTLIERS', 'REC', 'ALG'], ax_rec, full_list_algs, colors, 1.0)

    fig_f1 = plt.figure('F1')
    ax_f1 = fig_f1.add_axes([0.1, 0.1, 0.8, 0.8])
    #daltonic_boxplot(df, ['OUTLIERS', 'F1', 'ALG'], ax_f1, full_list_algs, colors, markers, 1.005)
    classic_boxplot(df, ['OUTLIERS', 'F1', 'ALG'], ax_f1, full_list_algs, colors, 1.0)

    fig_t = plt.figure('TIME')
    ax_t = fig_t.add_axes([0.1, 0.1, 0.8, 0.8])
    #daltonic_boxplot(df, ['OUTLIERS', 'TIME', 'ALG'], ax_ate, full_list_algs, colors, markers, -5.0, turn_off=True)
    classic_boxplot(df, ['OUTLIERS', 'TIME', 'ALG'], ax_t, full_list_algs, colors, -5.0)
    plt.ylim(0.01, 1000.0)
    plt.yscale('log')

    fig_ate = plt.figure('ATE')
    ax_ate = fig_ate.add_axes([0.1, 0.1, 0.8, 0.8])
    classic_boxplot(df, ['OUTLIERS', 'ATE', 'ALG'], ax_ate, full_list_algs, colors, 1.01)
    plt.ylim(0.1, 1000.0)
    plt.yscale('log')

    fig_rpe = plt.figure('RPE')
    ax_rpe = fig_rpe.add_axes([0.1, 0.1, 0.8, 0.8])
    #daltonic_boxplot(df, ['OUTLIERS', 'RPE', 'ALG'], ax_rpe, full_list_algs, colors, markers, 1.01, turn_off=True)
    classic_boxplot(df, ['OUTLIERS', 'RPE', 'ALG'], ax_rpe, full_list_algs, colors, 1.01)
    plt.ylim(0.01, 1000.0)
    plt.yscale('log')
    '''

    plt.show()

    return 



if __name__ == '__main__' :
    main()