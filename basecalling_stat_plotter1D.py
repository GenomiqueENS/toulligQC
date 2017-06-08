import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import re
import getter1D
import os



class basecalling_stat_plotter1D:
    """
    Plots different graphs for exploitation of minion runs from Albacore file log
    """

    def __init__(self, path_sequencing_summary, pdf, run_name, barcode_present, file_list=''):
        self.albacore_log = pd.read_csv(path_sequencing_summary, sep="\t")
        self.channel = self.albacore_log['channel']
        self.albacore_log[self.albacore_log == 0] = np.nan

        if barcode_present == 'y':
            self.barcode_selection = getter1D.get_Barcodes()
            self.barcode_selection_original = self.barcode_selection
            self.fast5_tot = len(self.albacore_log)
            self.pdf = pdf
            self.barcode_selection.sort()
            self.fastq_length_array = getter1D.get_FastqSeq_barcoded(self.barcode_selection, run_name, barcode_present, file_list)

        else:
            self.fast5_tot = len(self.albacore_log)
            self.pdf = pdf
            self.template_nucleotide_counter, self.total_nucs_template, self.fastq_length_array = getter1D.get_FastqSeq_without_barcode(run_name, file_list)

    def barcode_meanqscore(self):
        """
        Writes the mean qscore extracted from the log file provided by albacore
        """
        barcode_meanqscore_dataframe = self.albacore_log[['mean_qscore_template','barcode_arrangement']]
        barcode_selection = barcode_meanqscore_dataframe[barcode_meanqscore_dataframe['barcode_arrangement'].isin(self.barcode_selection)]

        for barcode in self.barcode_selection:
            meanq_score = barcode_selection[barcode_selection['barcode_arrangement'] == barcode].mean()
            meanq_score = meanq_score.tolist()[0]
            completeName = os.path.join('statistics/', barcode)
            file = open(completeName,'a')
            file.write("mean.qscore.template={}\n".format(meanq_score))
            file.close()
            print('fini')


    def run_date(self):
        """
        Returns the date of a Minion run from the log file provided by albacore
        """
        filename = self.albacore_log['filename']
        for index, file in enumerate(filename):
            exp = self.albacore_log['filename'][index]
            break
        m = re.search(r'(_(\d+)_)', exp)
        return m.group(2)

    def stat_generation(self):
        """
        Generates a dictionary of statistics such as quartile, the standard deviation for the creation of a log file from the log file provided by Albacore
        """
        num_called_template = self.albacore_log['num_called_template']
        mean_qscore_template = self.albacore_log['mean_qscore_template']
        statistics_num_called_template = pd.DataFrame.describe(num_called_template).drop("count")
        statistics_mean_qscore_template = pd.DataFrame.describe(mean_qscore_template).drop("count")
        return statistics_num_called_template, statistics_mean_qscore_template

    def barcode_percentage_pie_chart(self):
        """
        Plots a pie chart of the barcode percentage of a run. Needs the design file describing the barcodes to run
        """
        # Ne doit pas excéder 10

        for element in self.barcode_selection:


            if all(self.albacore_log['barcode_arrangement'] != element):
                print("The barcode {} doesn't exist".format(element))
                return False

        barcode = self.albacore_log['barcode_arrangement']
        count1 = barcode.value_counts()
        count = count1.sort_index()[self.barcode_selection]
        unclassified = sum(count1[~count1.index.isin(self.barcode_selection)])
        ##Watch out must be placed after the operation
        self.barcode_selection.append("unclassified")
        ##
        count['unclassified'] = unclassified
        total = sum(count1)

        cs = cm.Paired(np.arange(len(self.barcode_selection)) / len(self.barcode_selection))

        sizes = [(100 * chiffre) / total for chiffre in count.values]
        if len(self.barcode_selection) <= 10:
            fig1, ax1 = plt.subplots()
            ax1.pie(sizes, labels=self.barcode_selection, autopct='%1.1f%%', startangle=90, colors=cs)
            ax1.axis('equal')

        else:
            fig = plt.figure(figsize=(20, 10))
            ax1 = fig.add_subplot(111)
            length = np.arange(0, len(count))
            ax1.bar(length, count, color=cs)
            ax1.set_xticks(length)
            ax1.set_xticklabels(self.barcode_selection)
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image5.png')
        else:
            plt.savefig('images/image5.png')
        self.pdf.savefig()
        plt.close()


    #Launch after barcode pie chart because of self.barcode_selection
    def read_length_histogram(self):
        """
        Plots an histogram of the reads length by bins of 100 for each of the barcodes described in the design file or without barcode
        """
        if self.fastq_length_array == []:
            print('There is a mistake')
            return False
        plt.hist(self.fastq_length_array, edgecolor="#E6E6E6", color="#EE6666", bins=range(min(self.fastq_length_array), max(self.fastq_length_array) + 100, 100))
        plt.xlabel("fastq size(pb)")
        plt.ylabel("Count")
        plt.xlim(0,6000)
        plt.title("read size")
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image7.png')
        else:
            plt.savefig('images/image7.png')
        self.pdf.savefig()
        plt.close()

    def counter(self):
        """
        Paricipates to the count df nucleotide(A,T,C,G)" \
        """
        return self.template_nucleotide_counter, self.total_nucs_template

    def statistics_read_size(self):
        """
        Gets statistics from the file containing the fastq bz2 files decompressed
        """
        series_read_size = pd.Series(self.fastq_length_array)
        statistics = pd.Series.describe(series_read_size)
        return statistics


    def read_count_histogram(self):
        """
        Plots the count histograms of count  of the different types of reads eventually available in a Minion run: template, complement, full_2D.
        """

        # Count of fast5 total
        fast5tot = len(self.albacore_log)
        # Count of template reads
        template = len(self.albacore_log['num_called_template'].dropna())

         #Count of complement reads
        read_type = [fast5tot, template]
        label = ("fast5tot", "template")
        nd = np.arange(len(read_type))

        # Histogram of differents reads count(template, complement, fast2D)
        plt.bar(nd, read_type, align='center', color=["lightblue", "salmon"])
        plt.xticks(nd, label)
        plt.xlabel("read type")
        plt.ylabel("Counts")
        plt.title("Counts of read template")
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image1.png')
        else:
            plt.savefig('images/image1.png')
        self.pdf.savefig()
        plt.close()

    def read_quality_boxplot(self):
        """
        Plots a boxplot of reads quality
        """
        dataframe = self.albacore_log.loc[:, ["mean_qscore_template"]]
        sns.boxplot(data=dataframe)
        plt.title('Boxplot of read quality')
        plt.ylabel('Phred score')
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image2.png')
        else:
            plt.savefig('images/image2.png')
        self.pdf.savefig()
        plt.close()



    def channel_count_histogram(self):
        """
        Plots an histogram of the channel count according to the channel number
        """
        fig, ax = plt.subplots()
        ax.hist(self.albacore_log['channel'], edgecolor='black',  bins=range(min(self.albacore_log['channel']), max(self.albacore_log['channel']) + 64, 64))
        ax.set_xlabel("Channel number")
        ax.set_ylabel("Count")
        ax.set_title("Channel counts")
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image3.png')
        else:
            plt.savefig('images/image3.png')
        self.pdf.savefig()
        plt.close()



    def read_number_run(self):
        """
        Plots the reads produced along the run against the time(in hour)
        """
        start_time = self.albacore_log["start_time"] / 3600
        start_time_sorted = sorted(start_time)
        plt.scatter(start_time_sorted, np.arange(len(start_time_sorted)))
        plt.ylabel("produced reads")
        plt.xlabel("hour")
        plt.title("Read produced along the run")
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image4.png')
        else:
            plt.savefig('images/image4.png')
        self.pdf.savefig()
        plt.close()

    def minion_flowcell_layout(self):
        """
        Represents the layout of a minion flowcell
        """
        seeds = [125, 121, 117, 113, 109, 105, 101, 97,
                 93, 89, 85, 81, 77, 73, 69, 65,
                 61, 57, 53, 49, 45, 41, 37, 33,
                 29, 25, 21, 17, 13, 9, 5, 1]

        flowcell_layout = []
        for s in seeds:
            for block in range(4):
                 for row in range(4):
                    flowcell_layout.append(s + 128 * block + row)
        return flowcell_layout

    def plot_performance(self, pore_measure):
        """
        Plots the channels occupancy by the reads
        @:param pore_measure: reads number per pore
        """
        flowcell_layout = self.minion_flowcell_layout()

        pore_values = []
        for pore in flowcell_layout:
            if pore in pore_measure:
                pore_values.append(pore_measure[pore])
            else:
                pore_values.append(0)

        # make a data frame of the lists
        d = {'rownum': list(range(1, 17)) * 32,
             'colnum': sorted(list(range(1, 33)) * 16),
             'tot_reads': pore_values,
             'labels': flowcell_layout}

        df = pd.DataFrame(d)

        d = df.pivot("rownum", "colnum", "tot_reads")
        plt.figure(figsize=(20, 10))
        sns.heatmap(d, annot=True, fmt="d", linewidths=.5, cmap="YlGnBu")
        plt.title('Channel occupancy')
        if os.path.isfile('/configpass/docker_config.txt'):
            plt.savefig('/working.directory/images/image6.png')
        else:
            plt.savefig('images/image6.png')
        self.pdf.savefig()
        plt.close()


    def occupancy_pore(self):
        channel_count = self.channel
        total_number_reads_per_pore = pd.value_counts(channel_count)
        Series = pd.DataFrame.describe(total_number_reads_per_pore)
        return pd.Series.to_dict(Series)

    def get_barcode_selection(self):
        """
        Returns the selection of barcodes used from the design file.
        """
        return self.barcode_selection

    def statistics_dataframe(self):
        """
        Returns the statistics retrieved from the statistics files in the statistics directory for each barcode as a dataframe to make
        the reading easier.
        """
        df = pd.DataFrame(columns=self.barcode_selection)
        for selected_barcode in self.barcode_selection[:-1]:
            dico = {}
            file = open('statistics/{}'.format(selected_barcode),'r')
            for line in file:
                key, value = line.strip().split('=')
                dico[key.strip()] = value.strip()
            file.close()
            df[selected_barcode] = pd.Series(dico)
        df.to_csv('/home/ferrato/ownCloud/fast5_1D/dataframe.csv', header=self.selection_original,index=list(df.index), sep='\t')

