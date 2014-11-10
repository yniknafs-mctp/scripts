import numpy as np
import matplotlib.pyplot as plt


a = [0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,33,3,3,4,4,4,45,5,9]
b = [8,8,8,8,8,9,9,9,9,6,6,4,4,43,2, 200, 0 ,0]
c = [12,13,2,4,4,5,6,7,7,4,9, 100, 0 ,0]



def plotter(cnv, null, name, output_file, log_in = False):
    big_list = []
    
    for x in cnv:
        big_list.append([x, 'a', 'r', 1])
    for x in null:
        big_list.append([x, 'b', 'b', 0])
    
    #for x in normal:
    #    big_list.append([x, 'c', 'b'])
    big_list.sort(key=lambda k: (k[0]))
    
    expression_list = []
    progression_list = []
    color_list = []
    tick_list = []
    
    i = 0
    for x in big_list:
        i+=1
        expression_list.append(x[0])
        progression_list.append(x[1])
        color_list.append(x[2])
        if x[3] == 1:
            tick_list.append(i)
    ind = np.arange(len(big_list))
    
    width = 1       # the width of the bars: can also be len(x) sequence
    MIN_CONSTANT = 1+1e-3
    new_expression_list = np.log2(np.array(expression_list) + MIN_CONSTANT)
    plt.figure(figsize=(10,5))
    p1 = plt.bar(ind, new_expression_list, width, color=color_list, edgecolor='none')
    p2 = plt.bar([0], [MIN_CONSTANT], color='r', edgecolor='none')
    #p3 = plt.bar([0], [MIN_CONSTANT], color='y', edgecolor='none', log=log_in)
    p3 = plt.bar([0], [MIN_CONSTANT], color='b', edgecolor='none')
    
    def label_maker(label):
        return label + (width/2)
     
    plt.title(name)
    plt.ylabel('FPKM')
    plt.tick_params(labelsize='x-small')
    plt.xticks([])    
    #plt.yticks(np.arange(0,(1.25*max(expression_list)),(max(expression_list))/5))
    plt.xlabel('Total Samples (N=%d), CNV_samples (N=%d)' % (len(expression_list), len(cnv)))
    plt.xlim((0,1.01*len(expression_list)))
    plt.ylim((0,1.05*max(new_expression_list)))
    #plt.show()
    plt.vlines(tick_list, 0, 1)
    plt.legend((p2[0], p3[0]), ('CNV', 'Null'), loc=2, prop={'size':'small'})
    #plt.savefig(output_file)
    plt.show()
    
    #plt.yscale('log')

if __name__ == '__main__':
    pos = open('test_inputs_plot/pos.txt', 'r')
    for line in pos:
        pos_list = map(float, line.strip().replace(' ','').replace('[','').replace(']','').split(','))
    null = open('test_inputs_plot/null.txt', 'r')
    for line in null:
        null_list = map(float, line.strip().replace(' ','').replace('[','').replace(']','').split(','))
    norm = open('test_inputs_plot/norm.txt', 'r')
    for line in norm:
        norm_list = map(float, line.strip().replace(' ','').replace('[','').replace(']','').split(','))
         
    plotter(pos_list, null_list, norm_list, 'test', 'tester.png', True)
    new_list = pos_list + null_list + norm_list
    new_list_filter = filter(lambda a: a != 0, new_list)
    print min(new_list_filter)