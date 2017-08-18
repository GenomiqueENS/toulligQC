import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import re
import os
from toulligqc import fastq
import csv
from matplotlib import gridspec

from pandas.tools.plotting import table

class basecalling_stat_plotter1D:
    """
    Plots different graphs for exploitation of minion runs from Albacore file log
    """

    def __init__(self, path_sequencing_summary, barcode_present, result_directory, fastq_directory, fast5_tuple, run_name, is_barcode, fastq_file_extension, barcode_selection = '' ):


        self.minknown_version, self.flowcell_id, self.hostname, self.numMinion, self.run_id = fast5_tuple
        self.global_dictionnary = {}
        self.statistics_dictionnary = {}
        self.albacore_log = pd.read_csv(path_sequencing_summary, sep="\t")
        self.result_directory = result_directory
        self.channel = self.albacore_log['channel']
        self.sequence_length_template = self.albacore_log['sequence_length_template']
        self.null_event = self.albacore_log[self.albacore_log['num_events']==0]
        self.albacore_log = self.albacore_log.replace([np.inf, -np.inf], 0)
        self.albacore_log = self.albacore_log[self.albacore_log['num_events']!=0]
        fastq_object = fastq.fastq(result_directory, fastq_directory, run_name, is_barcode, fastq_file_extension)
        self.fast5_tot = len(self.albacore_log)
        if barcode_present:

            self.barcode_selection = barcode_selection
            self.albacore_log.loc[~self.albacore_log['barcode_arrangement'].isin(self.barcode_selection), 'barcode_arrangement'] = 'unclassified'
            self.barcode_selection.append('unclassified')
            self.fastq_length_array, self.global_dictionnary = fastq_object.get_fastq_barcoded(self.barcode_selection)
        else:
            self.total_nucs_template, self.fastq_length_array,_ , self.template_nucleotide_counter \
                = fastq_object.get_fastq_without_barcode()


    def make_table(self, value, ax, metric_suppression = ''):
        '''
        Creation of a statistics table printed with the graph
        :param value: information measured
        :param ax: axes used
        :param metric_suppression: suppression of a metric when we use the describe pandas function
        '''

        if metric_suppression:
            the_table = table(ax, np.round(value.describe().drop(metric_suppression), 2), loc='center', )
        else:
            the_table = table(ax, np.round(value.describe(), 2), loc='center', )

        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.axis('off')

        the_table.set_fontsize(12)
        the_table.scale(1, 1.2)

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
            file.write("mean.phred_score.template={}\n".format(meanq_score))
            file.close()


    def run_date(self):
        """
        Returns the date of a Minion run from the log file provided by albacore
        """
        file_name = self.albacore_log['filename'].iloc[0]
        pattern = re.search(r'(_(\d+)_)', file_name)
        return pattern.group(2)

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

        for element in self.barcode_selection:

            if all(self.albacore_log['barcode_arrangement'] != element):
                print("The barcode {} doesn't exist".format(element))
                return False

        barcode = self.albacore_log['barcode_arrangement']
        barcode_count = barcode.value_counts()
        count_sorted = barcode_count.sort_index()[self.barcode_selection]
        total = sum(count_sorted)

        cs = cm.Paired(np.arange(len(self.barcode_selection)) / len(self.barcode_selection))

        sizes = [(100 * chiffre) / total for chiffre in count_sorted.values]
        if len(self.barcode_selection) <= 10:
            fig1, ax1 = plt.subplots()
            ax1.pie(sizes, labels=self.barcode_selection, autopct='%.4f', startangle=90, colors=cs)
            ax1.axis('equal')

        else:
            fig = plt.figure(figsize=(20, 10))
            ax1 = fig.add_subplot(111)
            length = np.arange(0, len(barcode_count))
            ax1.bar(length, barcode_count, color=cs)
            ax1.set_xticks(length)
            ax1.set_xticklabels(self.barcode_selection)

        plt.savefig(self.result_directory+'images/barcode_percentage_pie_chart.png')
        plt.close()


    def safe_log(self,x):
        '''
        Verification that we haven't a null value
        :param x: tested value
        :return: log2 value or 0
        '''
        if x <= 0:
            return 0
        return np.log2(x)

    def read_length_histogram(self):
        """
        Plots an histogram of the reads length by bins of 100 for each of the barcodes described in the design file or without barcode
        """

        minimum, maximum = min(self.sequence_length_template), max(self.sequence_length_template)

        fig = plt.figure(figsize=(20, 10))
        gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
        ax = plt.subplot(gs[0])

        n, bins, patches = ax.hist(self.sequence_length_template,edgecolor ='black', bins= 2**np.linspace(self.safe_log(minimum), self.safe_log(maximum),30))
        ax.set_xscale('log',basex=2)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
        ax.set_xticks(bins)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_rotation('vertical')
        ax.set_xlabel('Read size(in pb)')
        ax.set_ylabel('Read number')
        ax.set_title('Read size histogram')

        ax2 = plt.subplot(gs[1])

        the_table = table(ax2, np.round(self.sequence_length_template.describe(), 2), loc='center')

        ax2.xaxis.set_visible(False)
        ax2.yaxis.set_visible(False)
        ax2.axis('off')

        the_table.set_fontsize(12)
        the_table.scale(1, 2)

        plt.savefig(self.result_directory + 'images/read_length_histogram.png')
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

    def phred_score_frequency(self):
        '''
        Plot the distribution of the phred score
        '''
        figure = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])

        sns.distplot(self.albacore_log['mean_qscore_template'], bins=15, color='green',
                     hist_kws=dict(edgecolor="k", linewidth=1))
        plt.xlabel("mean_qscore")
        plt.ylabel("Frequency")
        plt.title("Phred score frequency")
        rd = self.albacore_log['mean_qscore_template'].describe().drop('count').round(2).reset_index()

        plt.axvline(x=self.albacore_log['mean_qscore_template'].describe()['50%'], color='green')

        ax2 = plt.subplot(gs[1])
        self.make_table(rd, ax2, 'count')
        plt.savefig(self.result_directory + 'images/phred_score_frequency.png')
        plt.close()

    def scatterplot(self):
        '''
        Plot the scatter plot representing the relation between the phred score and the sequence length
        '''
        plt.scatter(x = self.albacore_log['sequence_length_template'], y = self.albacore_log['mean_qscore_template'])
        plt.xlim(0,100000)
        plt.xlabel("sequence_length_template")
        plt.ylabel("mean_qscore_template")
        plt.title("Relation between the sequence length template and the mean qscore template")
        plt.savefig(self.result_directory + 'images/scatter_plot.png')
        plt.close()

    def read_count_histogram(self):
        """
        Plots the count histograms of count  of the different types of reads eventually available in a Minion run: template, complement, full_2D.
        """

        fast5_raw = len(self.albacore_log['num_events'])
        fast5_template_basecalled = \
            len(self.albacore_log['num_events'])-len(self.albacore_log[self.albacore_log['num_called_template']==0])
        self.statistics_dictionnary['fast5_template_NOT_basecalled'] = \
            self.albacore_log[self.albacore_log['num_called_template'] == 0]
        self.statistics_dictionnary['fast5_template_NOT_basecalled_percentage'] = \
            (len(self.albacore_log[self.albacore_log['num_called_template']==0]))/len(self.albacore_log['num_called_template'])*100

        read_type = [fast5_raw, fast5_template_basecalled]
        label = ("fast5_raw", "fast5_template_basecalled")
        nd = np.arange(len(read_type))

        bars = plt.bar(nd, read_type, align='center', color=["lightblue", "salmon"])
        plt.xticks(nd, label)
        plt.xlabel("read type")
        plt.ylabel("Counts")
        plt.title("Counts of read template")


        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2., 1 * height, '%d' % int(height), ha='center', va='bottom')

        plt.savefig(self.result_directory+'images/read_count_histogram.png')
        plt.close()

    def read_quality_boxplot(self):
        """
        Plots a boxplot of reads quality
        """
        figure = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])

        dataframe = self.albacore_log.loc[:, ["mean_qscore_template"]]
        sns.boxplot(data=dataframe)
        plt.title('Boxplot of read quality')
        plt.ylabel('Phred score')

        ax2 = plt.subplot(gs[1])
        self.make_table(dataframe, ax2, 'count')

        plt.savefig(self.result_directory+'images/read_quality_boxplot.png')
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

        plt.savefig(self.result_directory+'images/channel_count_histogram.png')
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

        plt.savefig(self.result_directory+'images/read_number_run.png')
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

        d = {'rownum': list(range(1, 17)) * 32,
             'colnum': sorted(list(range(1, 33)) * 16),
             'tot_reads': pore_values,
             'labels': flowcell_layout}

        df = pd.DataFrame(d)

        d = df.pivot("rownum", "colnum", "tot_reads")
        d2 = df.pivot("rownum", "colnum","labels")
        plt.figure(figsize=(20, 10))
        sns.heatmap(d, fmt="", annot = d2, linewidths=.5, cmap="YlGnBu")
        plt.title('Channel occupancy')

        plt.savefig(self.result_directory+'images/channel_occupancy.png')
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


    def statistics_dictionnary(self):
        nucleotide = ['A', 'T', 'C', 'G']
        sequence_length_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
        channel_occupancy_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
        mean_qscore_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
        nucleotide_count_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(4)]
        nucleotide_proportion_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(4)]


        column = 0
        for index_barcode, barcode in enumerate(self.barcode_selection):
            barcode_selected_dataframe = self.albacore_log[self.albacore_log['barcode_arrangement'] == barcode]
            channel_occupancy_statistics = barcode_selected_dataframe['channel'].describe()
            mean_qscore_statistics = barcode_selected_dataframe['mean_qscore_template'].describe()
            sequence_length_statistics = barcode_selected_dataframe['sequence_length_template'].describe()
            sorted_list = sorted(self.global_dictionnary['nucleotide_count_'+barcode].items())
            if column == 0:
                for i in range(8):
                    mean_qscore_matrix[i][column] = 'phred_score_' + mean_qscore_statistics.keys()[i]
                    channel_occupancy_matrix[i][0] = 'channel_occupancy_' + channel_occupancy_statistics.keys()[i]
                    sequence_length_matrix[i][0] = 'sequence_length_' + sequence_length_statistics.keys()[i]

                for line, nucleotide_number_list in enumerate(sorted_list):
                    if nucleotide_number_list[0] in nucleotide:
                        nucleotide_count_matrix[line][column] = 'nucleotide_count_' + nucleotide_number_list[0]
                        nucleotide_proportion_matrix[line][column] = 'nucleotide_proportion_'+ nucleotide_number_list[0]
                    else:
                        continue

                column += 1

            for line, metric in enumerate(mean_qscore_statistics):
                mean_qscore_matrix[line][column] = round(metric, 3)

            for line, metric in enumerate(channel_occupancy_statistics):
                channel_occupancy_matrix[line][column] = round(metric, 3)

            for line, metric in enumerate(sequence_length_statistics):
                sequence_length_matrix[line][column] = round(metric, 3)

            for line, nucleotide_number_list in enumerate(sorted_list):
                if nucleotide_number_list[0] in nucleotide:
                    nucleotide_count_matrix[line][column] = nucleotide_number_list[1]
                    calcul = nucleotide_number_list[1] / self.global_dictionnary['total_nucleotide_'+barcode]
                    nucleotide_proportion_matrix[line][column] = calcul
                else:
                    continue
            column += 1

    def statistics_dataframe(self):
        """
        Returns the statistics retrieved from the statistics files in the statistics directory for each barcode as a dataframe to make
        the reading easier.
        """

        nucleotide = ['A', 'T', 'C', 'G']
        sequence_length_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
        channel_occupancy_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
        mean_qscore_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(8)]
        nucleotide_count_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(4)]
        nucleotide_proportion_matrix = [[0 for column in range(len(self.barcode_selection) + 1)] for line in range(4)]


        column = 0
        for index_barcode, barcode in enumerate(self.barcode_selection):
            barcode_selected_dataframe = self.albacore_log[self.albacore_log['barcode_arrangement'] == barcode]
            channel_occupancy_statistics = barcode_selected_dataframe['channel'].describe()
            mean_qscore_statistics = barcode_selected_dataframe['mean_qscore_template'].describe()
            sequence_length_statistics = barcode_selected_dataframe['sequence_length_template'].describe()
            if barcode != 'unclassified':
                sorted_list = sorted(self.global_dictionnary['nucleotide_count_'+barcode].items())
            if column == 0:
                for i in range(8):
                    mean_qscore_matrix[i][column] = 'phred_score_' + mean_qscore_statistics.keys()[i]
                    channel_occupancy_matrix[i][0] = 'channel_occupancy_' + channel_occupancy_statistics.keys()[i]
                    sequence_length_matrix[i][0] = 'sequence_length_' + sequence_length_statistics.keys()[i]

                for line, nucleotide_number_list in enumerate(sorted_list):
                    if nucleotide_number_list[0] in nucleotide and barcode != 'unclassified':
                        nucleotide_count_matrix[line][column] = 'nucleotide_count_' + nucleotide_number_list[0]
                        nucleotide_proportion_matrix[line][column] = 'nucleotide_proportion_'+ nucleotide_number_list[0]
                    else:
                        continue

                column += 1

            for line, metric in enumerate(mean_qscore_statistics):
                mean_qscore_matrix[line][column] = round(metric, 3)

            for line, metric in enumerate(channel_occupancy_statistics):
                channel_occupancy_matrix[line][column] = round(metric, 3)

            for line, metric in enumerate(sequence_length_statistics):
                sequence_length_matrix[line][column] = round(metric, 3)

            for line, nucleotide_number_list in enumerate(sorted_list):
                if nucleotide_number_list[0] in nucleotide and barcode != 'unclassified':
                    nucleotide_count_matrix[line][column] = float(nucleotide_number_list[1])
                    calcul = nucleotide_number_list[1] / self.global_dictionnary['total_nucleotide_'+barcode]
                    nucleotide_proportion_matrix[line][column] = round(calcul, 3)
                else:
                    continue
            column += 1

        barcode_selection_matrix = list(self.barcode_selection)
        barcode_selection_matrix.insert(0, '')
        self.global_dictionnary["minknown.version"] = self.minknown_version
        self.global_dictionnary["hostname"] = self.hostname
        self.global_dictionnary["minion.serial.number"] = self.numMinion
        self.global_dictionnary["run.id"] = self.run_id
        general_information_list = [['minknown', self.global_dictionnary["minknown.version"]],
                                    ['hostname', self.global_dictionnary["hostname"]],
                                    ['minion.serial.number', self.global_dictionnary["minion.serial.number"]],
                                    ['run.id', self.global_dictionnary["run.id"]]]

        with open(self.result_directory + 'dataframe.csv', 'a') as csv_file:
            writer = csv.writer(csv_file, delimiter='\t')

            writer.writerow(['Number of reads:', len(self.sequence_length_template)])
            writer.writerow('')

            writer.writerow(barcode_selection_matrix)
            for element in channel_occupancy_matrix:
                writer.writerow(element)

            writer.writerow('')
            for metric in mean_qscore_matrix:
                writer.writerow(metric)

            writer.writerow('')
            for metric in sequence_length_matrix:
                writer.writerow(metric)

            writer.writerow('')
            for metric in nucleotide_count_matrix:
                writer.writerow(metric)

            writer.writerow('')
            for metric in nucleotide_proportion_matrix:
                writer.writerow(metric)

            writer.writerow('')
            writer.writerows(general_information_list)

    def barcode_length_boxplot(self):
        '''
        Plot the length boxplot for each barcode indicated in the sample sheet
        '''
        pattern = '(\d{2})'
        dico = {}


        fig = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])


        for barcode in self.barcode_selection:
            barcode_selected_dataframe = self.albacore_log[self.albacore_log['barcode_arrangement'] == barcode]
            match = re.search(pattern, barcode)
            if match:
                dico[match.group(0)] = barcode_selected_dataframe['sequence_length_template']
            else:
                dico['U'] = barcode_selected_dataframe['sequence_length_template']
        barcode_selection_sequence_length_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
        sns.boxplot(data=barcode_selection_sequence_length_dataframe, showfliers=False)
        plt.xlabel('Barcodes')
        plt.ylabel('Read size(in pb)')
        plt.title('Read size distribution for each barcode')

        ax2 = plt.subplot(gs[1])
        self.make_table(barcode_selection_sequence_length_dataframe,ax2)
        plt.savefig(self.result_directory + 'images/barcode_length_boxplot.png')
        plt.close()

    def barcoded_phred_score_frequency(self):
        '''
        Plot the phred score distribution boxplot for each barcode indicated in the sample sheet
        '''
        dico = {}
        pattern = '(\d{2})'

        fig = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])

        for barcode in self.barcode_selection:
            barcode_selected_dataframe = self.albacore_log[self.albacore_log['barcode_arrangement'] == barcode]
            match = re.search(pattern, barcode)
            if match:
                dico[match.group(0)] = barcode_selected_dataframe['mean_qscore_template']
            else:
                dico['U'] = barcode_selected_dataframe['mean_qscore_template']
        barcode_selection_phred_scrore_dataframe = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in dico.items()]))
        sns.boxplot(data= barcode_selection_phred_scrore_dataframe, showfliers=False)
        plt.xlabel('Barcodes')
        plt.ylabel('Phred score')
        plt.title('Phred score distribution for each barcode')

        ax2 = plt.subplot(gs[1])
        self.make_table(barcode_selection_phred_scrore_dataframe,ax2)
        plt.savefig(self.result_directory + 'images/barcode_phred_score_boxplot.png')
        plt.close()


