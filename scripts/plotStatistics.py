import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

def getMeanFromDict(dict_stats) :

    full_mean = np.array([])
    for elem in dict_stats :
        stats = dict_stats[elem]
        stats = np.reshape(stats, (3, -1))
        mean = np.mean(stats, axis=1)
        full_mean = np.append(full_mean, mean)

    return np.reshape(full_mean, (-1, 3))

def main() :

    # Reading the data
    with open('/home/slam-emix/Workspace/BackendBench/src/voting_consensus/scripts/data/M3500_20.txt') as f:
        lines = f.readlines()

    # EST | GT | N_VARIABELS | N_LOOPS | ITERS | DELTA_T

    vect = []
    colors_gt = []
    colors_est = []
    gt_off = []
    est_off = []
    diz1_n_variables = {}
    diz1_n_variables_out = {}
    diz1_n_variables_in = {}
    c = 0
    for line in lines :
        x = line.split()
        v = [ int(x[0]), int(x[1]), int(x[2]), 
             int(x[3]), int(x[4]), float(x[5][:-2])]
        vect.append(v)
        cl_gt = cl_est = "blue"
        val_gt = val_est = 0 
        if vect[-1][1] == 0 :
            cl_gt = "red"
            val_gt = -1
        if vect[-1][0] == 0 :
            cl_est = "red"
            val_est = -1

        colors_gt.append(cl_gt)
        colors_est.append(cl_est)
        gt_off.append(val_gt)
        est_off.append(val_est)

        if ( v[2] in diz1_n_variables ) :
            diz1_n_variables[v[2]] = np.append(diz1_n_variables[v[2]], [v[-3], v[-2], v[-1]])
            c += 1
        else :
            diz1_n_variables[v[2]] = np.array([v[-3], v[-2], v[-1]])

        if ( vect[-1][1] == 0 and v[2] in diz1_n_variables_out ) :
            diz1_n_variables_out[v[2]] = np.append(diz1_n_variables_out[v[2]], [v[-3], v[-2], v[-1]])
        elif ( vect[-1][1] == 0 and v[2] not in diz1_n_variables_out ) :
            diz1_n_variables_out[v[2]] = np.array([v[-3], v[-2], v[-1]])

        if ( vect[-1][1] == 1 and v[2] in diz1_n_variables_in ) :
            diz1_n_variables_in[v[2]] = np.append(diz1_n_variables_in[v[2]], [v[-3], v[-2], v[-1]])
        elif ( vect[-1][1] == 1 and v[2] not in diz1_n_variables_in ) :
            diz1_n_variables_in[v[2]] = np.array([v[-3], v[-2], v[-1]])
        
    vect = np.array(vect)
    time_d = np.arange(0, vect.shape[0], 1, dtype=int)
    dis_gt  = vect[:, 1] + np.array(gt_off)
    dis_est = vect[:, 0] + np.array(est_off)

    '''
    plt.figure(1)
    plt.subplot(211)
    plt.title("Probability Distribution")
    plt.ylabel("Distribution")
    plt.bar(time_d, np.ones(time_d.shape, dtype=int), color = colors_gt)
    plt.subplot(212)
    plt.bar(time_d, np.ones(time_d.shape, dtype=int), color = colors_est)
    plt.ylabel("Distribution")
    plt.xlabel("Time")

    cumulative_gt = []
    cumulative_est = []
    sum_gt = 0
    sum_est = 0
    for elem in vect :
        sum_est += elem[0]
        sum_gt += elem[1]
        cumulative_gt.append(sum_gt)
        cumulative_est.append(sum_est)

    cumulative_gt = np.array(cumulative_gt)
    cumulative_est = np.array(cumulative_est)
    plt.figure(2)
    plt.title("Cumulative Distribution")
    plt.xlabel("Time")
    plt.ylabel("Cumulative")
    plt.step(time_d, cumulative_gt, label='GT')
    plt.step(time_d, cumulative_est, label='EST')
    #plt.plot(time_d, cumulative_gt, label='GT')
    #plt.plot(time_d, cumulative_est, label='EST')
    plt.legend()
    

    plt.figure(3)
    plt.title("BoxPlots")
    plt.boxplot(np.stack((dis_est, dis_gt), axis=1), labels=['est', 'gt'])

    print('INFO EST')
    print("- mean: "+str(statistics.mean(dis_est)))
    print("- Q1: "+str(np.percentile(dis_est, 25)))
    print("- median: "+str(np.percentile(dis_est, 50)))
    print("- Q3: "+str(np.percentile(dis_est, 75)))

    print('INFO GT')
    print("- mean: "+str(statistics.mean(dis_gt)))
    print("- Q1: "+str(np.percentile(dis_gt, 25)))
    print("- median: "+str(np.percentile(dis_gt, 50)))
    print("- Q3: "+str(np.percentile(dis_gt, 75)))
    '''

    values_mixed = getMeanFromDict(diz1_n_variables)
    values_in  = getMeanFromDict(diz1_n_variables_in)
    values_out = getMeanFromDict(diz1_n_variables_out)
    
    plt.figure('Mixed')
    #plt.subplot(311)
    #plt.title("Relations over time")
    #plt.ylabel("Loops")
    #plt.bar(diz1_n_variables.keys(), values_mixed[:, 0])
    #plt.subplot(312)
    plt.ylabel("Iters")
    plt.bar(diz1_n_variables.keys(), values_mixed[:, 1])
    #plt.subplot(313)
    plt.xlabel("N Variables")
    #plt.ylabel("Time")
    #plt.bar(diz1_n_variables.keys(), values_mixed[:, 2])

    
    if (len(diz1_n_variables_in) > 0 ) :
        plt.figure('Inliers')
        plt.subplot(311)
        plt.title("Relations over time")
        plt.ylabel("Loops")
        plt.bar(diz1_n_variables_in.keys(), values_in[:, 0])
        plt.subplot(312)
        plt.ylabel("Iters")
        plt.bar(diz1_n_variables_in.keys(), values_in[:, 1])
        plt.subplot(313)
        plt.xlabel("N Variables")
        plt.ylabel("Time")
        plt.bar(diz1_n_variables_in.keys(), values_in[:, 2])

    if (len(diz1_n_variables_out) > 0 ) :
        plt.figure('Outliers')
        plt.subplot(311)
        plt.title("Relations over time")
        plt.ylabel("Loops")
        plt.bar(diz1_n_variables_out.keys(), values_out[:, 0])
        plt.subplot(312)
        plt.ylabel("Iters")
        plt.bar(diz1_n_variables_out.keys(), values_out[:, 1])
        plt.subplot(313)
        plt.xlabel("N Variables")
        plt.ylabel("Time")
        plt.bar(diz1_n_variables_out.keys(), values_out[:, 2])

    plt.figure('3D')
    ax = plt.axes(projection='3d')
    ax.scatter3D(list(diz1_n_variables.keys()), values_mixed[:, 0], values_mixed[:, 1])
    ax.set_xlabel('VARIABLES')
    ax.set_ylabel('LOOPS')
    ax.set_zlabel('ITERS')

    plt.show()

    return 



if __name__ == '__main__' :
    main()