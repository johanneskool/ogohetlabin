import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# Define a linear function for curve fitting
def linear(x, m, b):
    return m * x + b


def sigmoid(x, A, B, C, D):
    # 4-parameter sigmoid function
    return A + (B - A) / (1 + (x / C) ** D)


# Function to calculate the mean across lists for each row
def mean_of_lists(row):
    # Zip the lists together to align corresponding elements and calculate mean for each position
    zipped_lists = zip(*row)
    mean_list = [np.mean(elements) for elements in zipped_lists]
    return mean_list


def extract_sheet(xlsx):
    # Get the directory of the current script
    WORKING_DIR = os.path.dirname(os.path.abspath(__file__))
    # Define file paths
    file_path = os.path.join(WORKING_DIR, '%s' % xlsx)
    # load section into dataframe
    # headers are in row 44, labels in column A.
    # the data starts at row 45 and ends at row 142
    df = pd.read_excel(file_path, header=43, nrows=98, index_col=0, )
    # the data is now structured as follows:
    # rows are labeled A1-12, B1-12, ... H1-12
    # columns are labeled 1-12, and the data is in the cells as a number
    # Restructure the data such that rows are labeled A-H
    # and columns are labeled 1-12
    # the data must now be a sequence of 12 numbers
    concentrations = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    df2 = pd.DataFrame(index=concentrations, columns=columns)
    for row in concentrations:
        for column in columns:
            corresponding_row = df.loc[row + str(column)]
            df2.at[row, column] = corresponding_row.values.tolist()
    time_points = df.loc['Time [s]']
    temperatures = df.loc['Temp. [°C]']

    ## DIPEPTIDE 1
    # extract rows A-H in columns 4,5,6, but call them 1-3 in the new df
    df_dipep1_triplo = pd.DataFrame(index=concentrations, columns=[1, 2, 3])
    for row in concentrations:
        corresponding_row = df2.loc[row]
        df_dipep1_triplo.at[row, 1] = corresponding_row[4]
        df_dipep1_triplo.at[row, 2] = corresponding_row[5]
        df_dipep1_triplo.at[row, 3] = corresponding_row[6]
    df_dipep1_triplo['mean'] = df_dipep1_triplo.apply(mean_of_lists, axis=1)
    # create a new datafram with concentration as the index and the values at each timepoint as the columns
    df_dipep1 = pd.DataFrame(df_dipep1_triplo['mean'].tolist(), index=concentrations, columns=time_points)

    ## DIPEPTIDE 2
    df_dipep2_triplo = pd.DataFrame(index=concentrations, columns=[1, 2, 3])
    for row in concentrations:
        corresponding_row = df2.loc[row]
        df_dipep2_triplo.at[row, 1] = corresponding_row[7]
        df_dipep2_triplo.at[row, 2] = corresponding_row[8]
        df_dipep2_triplo.at[row, 3] = corresponding_row[9]
    df_dipep2_triplo['mean'] = df_dipep2_triplo.apply(mean_of_lists, axis=1)
    df_dipep2 = pd.DataFrame(df_dipep2_triplo['mean'].tolist(), index=concentrations, columns=time_points)

    ## DIPEPTIDE 3
    df_dipep3_triplo = pd.DataFrame(index=concentrations, columns=[1, 2, 3])
    for row in concentrations:
        corresponding_row = df2.loc[row]
        df_dipep3_triplo.at[row, 1] = corresponding_row[10]
        df_dipep3_triplo.at[row, 2] = corresponding_row[11]
        df_dipep3_triplo.at[row, 3] = corresponding_row[12]
    df_dipep3_triplo['mean'] = df_dipep3_triplo.apply(mean_of_lists, axis=1)
    df_dipep3 = pd.DataFrame(df_dipep3_triplo['mean'].tolist(), index=concentrations, columns=time_points)

    ## Captopril
    df_Captopril_triplo = pd.DataFrame(index=['A'], columns=[1, 2, 3])
    corresponding_row = df2.loc['A']
    df_Captopril_triplo.at['A', 1] = corresponding_row[1]
    df_Captopril_triplo.at['A', 2] = corresponding_row[2]
    df_Captopril_triplo.at['A', 3] = corresponding_row[3]
    df_Captopril_triplo['mean'] = df_Captopril_triplo.apply(mean_of_lists, axis=1)
    df_Captopril = pd.DataFrame(df_Captopril_triplo['mean'].tolist(), index=['A'], columns=time_points)

    ## GlyGly
    df_GlyGly_triplo = pd.DataFrame(index=['A'], columns=[1, 2, 3])
    corresponding_row = df2.loc['B']
    df_GlyGly_triplo.at['A', 1] = corresponding_row[1]
    df_GlyGly_triplo.at['A', 2] = corresponding_row[2]
    df_GlyGly_triplo.at['A', 3] = corresponding_row[3]
    df_GlyGly_triplo['mean'] = df_GlyGly_triplo.apply(mean_of_lists, axis=1)
    df_GlyGly = pd.DataFrame(df_GlyGly_triplo['mean'].tolist(), index=['A'], columns=time_points)

    ## Buff_Ace_Sub
    df_Buff_Ace_Sub_triplo = pd.DataFrame(index=['A'], columns=[1, 2, 3])
    corresponding_row = df2.loc['C']
    df_Buff_Ace_Sub_triplo.at['A', 1] = corresponding_row[1]
    df_Buff_Ace_Sub_triplo.at['A', 2] = corresponding_row[2]
    df_Buff_Ace_Sub_triplo.at['A', 3] = corresponding_row[3]
    df_Buff_Ace_Sub_triplo['mean'] = df_Buff_Ace_Sub_triplo.apply(mean_of_lists, axis=1)
    df_Buff_Ace_Sub = pd.DataFrame(df_Buff_Ace_Sub_triplo['mean'].tolist(), index=['A'], columns=time_points)

    ## Buff_Sub
    df_Buff_Sub_triplo = pd.DataFrame(index=['A'], columns=[1, 2, 3])
    corresponding_row = df2.loc['D']
    df_Buff_Sub_triplo.at['A', 1] = corresponding_row[1]
    df_Buff_Sub_triplo.at['A', 2] = corresponding_row[2]
    df_Buff_Sub_triplo.at['A', 3] = corresponding_row[3]
    df_Buff_Sub_triplo['mean'] = df_Buff_Sub_triplo.apply(mean_of_lists, axis=1)
    df_Buff_Sub = pd.DataFrame(df_Buff_Sub_triplo['mean'].tolist(), index=['A'], columns=time_points)

    ## Buff
    df_Buff_triplo = pd.DataFrame(index=['A'], columns=[1, 2, 3])
    corresponding_row = df2.loc['E']
    df_Buff_triplo.at['A', 1] = corresponding_row[1]
    df_Buff_triplo.at['A', 2] = corresponding_row[2]
    df_Buff_triplo.at['A', 3] = corresponding_row[3]
    df_Buff_triplo['mean'] = df_Buff_triplo.apply(mean_of_lists, axis=1)
    df_Buff = pd.DataFrame(df_Buff_triplo['mean'].tolist(), index=['A'], columns=time_points)

    return concentrations, time_points, df_dipep1, df_dipep2, df_dipep3, df_Captopril, df_GlyGly, df_Buff_Ace_Sub, df_Buff_Sub, df_Buff


def plot_dipeptides_with_trendline(datasets, time_points, df_Captopril, df_GlyGly, df_Buff_Ace_Sub):
    # dataframes for dipeptide slopes at concentrations A-H
    df_dipep1_slopes = pd.DataFrame(index=datasets[0][1].index, columns=[0])
    df_dipep2_slopes = pd.DataFrame(index=datasets[1][1].index, columns=[0])
    df_dipep3_slopes = pd.DataFrame(index=datasets[2][1].index, columns=[0])

    # setup colormap
    colormap = plt.cm.inferno
    # Generate 8 colors from the colormap
    # np.linspace(start, stop, num) generates 'num' evenly spaced samples, calculated over the interval [start, stop].
    conc_colors = dict(zip(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], [colormap(i) for i in np.linspace(0, 1, 8)]))

    # Create a single figure and a grid of subplots
    fig, axes = plt.subplots(nrows=1, ncols=len(datasets), figsize=(24, 8), sharey=True)  # Adjust figsize as needed
    # If there's only one subplot, axes might not be an array. Ensure it's iterable.
    if len(datasets) == 1:
        axes = [axes]

    labels = {'A': '2187 µM', 'B': '729 µM', 'C': '243 µM', 'D': '81 µM', 'E': '27 µM', 'F': '9 µM', 'G': '3 µM',
              'H': '1 µM'}

    dipep_num = 1
    for ax, (name, df) in zip(axes, datasets):
        # plot MetLeu inhibition at concentrations A-H over time
        for row in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']:  # concentrations A-H
            # for row in ['A', 'H']:  # concentrations A-H
            label = "[{}]".format(labels[row])
            color = conc_colors[row]

            inhibition_values = df.loc[row]

            fit_linear(ax, color, inhibition_values, time_points)
            # add equation to label
            curve = curve_fit(linear, np.array(time_points), np.array(inhibition_values))
            if dipep_num == 1:
                df_dipep1_slopes.loc[row] = curve[0][0]
            elif dipep_num == 2:
                df_dipep2_slopes.loc[row] = curve[0][0]
            elif dipep_num == 3:
                df_dipep3_slopes.loc[row] = curve[0][0]

            label += ' y = {:.2f}x + {:.2f}'.format(*curve[0])

            ax.plot(time_points, inhibition_values, 'x', label=label, color=color)

        # compare against Captopril
        fit_linear(ax, 'red', df_Captopril.loc['A'], time_points)
        label_Captopril = 'Captopril 10 µM'

        curve_Captopril = curve_fit(linear, np.array(time_points), np.array(df_Captopril.loc['A']))
        label_Captopril += ' y = {:.2f}x + {:.2f}'.format(*curve_Captopril[0])
        Captopril_slope = curve_Captopril[0][0]
        ax.plot(time_points, df_Captopril.loc['A'], 'x', label=label_Captopril, color='red')

        # compare against GlyGly
        fit_linear(ax, 'blue', df_GlyGly.loc['A'], time_points)
        label_GlyGly = 'GlyGly 10 µM'
        curve_GlyGly = curve_fit(linear, np.array(time_points), np.array(df_GlyGly.loc['A']))
        label_GlyGly += ' y = {:.2f}x + {:.2f}'.format(*curve_GlyGly[0])
        GlyGly_slope = curve_GlyGly[0][0]
        ax.plot(time_points, df_GlyGly.loc['A'], 'x', label=label_GlyGly, color='blue')

        # compare against Buff_Ace_Sub
        fit_linear(ax, 'green', df_Buff_Ace_Sub.loc['A'], time_points)
        label_Buff_Ace_Sub = 'Buff_Ace_Sub'
        curve_Buff_Ace_Sub = curve_fit(linear, np.array(time_points), np.array(df_Buff_Ace_Sub.loc['A']))
        label_Buff_Ace_Sub += ' y = {:.2f}x + {:.2f}'.format(*curve_Buff_Ace_Sub[0])
        Buff_Ace_Sub_slope = curve_Buff_Ace_Sub[0][0]
        ax.plot(time_points, df_Buff_Ace_Sub.loc['A'], 'x', label=label_Buff_Ace_Sub, color='green')

        ax.set_xlabel('Time [s]')
        ax.set_xticks(time_points)
        ax.set_title(f'Inhibition over time for {name}')

        ax.legend()

        dipep_num += 1
    # Only add the y-label to the leftmost subplot for clarity
    axes[0].set_ylabel('inhibition')
    fig.savefig("{}-{}-{}.png".format(datasets[0][0], datasets[1][0], datasets[2][0]))
    # plt.tight_layout()
    plt.show()

    return df_dipep1_slopes, df_dipep2_slopes, df_dipep3_slopes, Captopril_slope, GlyGly_slope, Buff_Ace_Sub_slope


def fit_linear(ax, color, inhibition_values, time_points):
    # fit linear region with a trendline
    x = np.array(time_points)
    y = np.array(inhibition_values)
    params, _ = curve_fit(linear, x, y)
    # Use the fitted parameters to plot the trendline
    trendline = linear(x, *params)
    ax.plot(x, trendline, '--', color=color)  # Dashed line for the trendline


def fit_sigmoid(ax, color, concentrations, slopes):
    # fit a sigmoid function to the data
    x = np.log10(concentrations)  # Use log10 if plotting on a log scale    y = np.array(slopes)
    y = slopes
    try:
        # Provide initial parameter guess and set parameter bounds
        # Set initial guesses to reasonable values based on your data
        p0 = [min(y), max(y), np.median(x), 1]
        bounds = ([0, 0, -np.inf, 0], [np.inf, np.inf, np.inf, np.inf])  # Example bounds

        params, cov = curve_fit(sigmoid, x, y, p0=p0, bounds=bounds)

        # Generate x values for the sigmoid curve
        x_fit = np.linspace(min(x), max(x), 100)
        trendline = sigmoid(x_fit, *params)

        # Plot the trendline
        ax.plot(10 ** x_fit, trendline, '--', color=color)  # Convert back if x was logged
    except RuntimeError as e:
        print(f"Could not fit sigmoid to data: {e}")


if __name__ == '__main__':
    concentrations1, time_points1, df_LeuIle, df_LeuLeu, df_LeuVal, df_Captopril1, df_GlyGly1, df_Buff_Ace_Sub1, df_Buff_Sub1, df_Buff1 = extract_sheet(
        xlsx='data/raw/data-1-2-6.xlsx')

    concentrations2, time_points2, df_MetIle, df_MetLeu, df_MetVal, df_Captopril2, df_GlyGly2, df_Buff_Ace_Sub2, df_Buff_Sub2, df_Buff2 = extract_sheet(
        xlsx='data/raw/data-3-4-5.xlsx')

    ### PLOTTING ###

    dipeptides1 = [('LeuIle', df_LeuIle), ('LeuLeu', df_LeuLeu), ('LeuVal', df_LeuVal)]
    # plot_dipeptides_with_trendline(dipeptides1, time_points1, df_Captopril1, df_GlyGly1)
    df_LeuIle_slopes, df_LeuLeu_slopes, df_LeuVal_slopes, Captopril1_slope, GlyGly1_slope, Buff_Ace_Sub1_slope = plot_dipeptides_with_trendline(
        dipeptides1, time_points1, df_Captopril1, df_GlyGly1, df_Buff_Ace_Sub1)

    dipeptides2 = [('MetVal', df_MetVal), ('MetIle', df_MetIle), ('MetLeu', df_MetLeu)]
    # plot_dipeptides_with_trendline(dipeptides2, time_points2, df_Captopril2, df_GlyGly2)
    df_MetVal_slopes, df_MetIle_slopes, df_MetLeu_slopes, Captopril2_slope, GlyGly2_slope, Buff_Ace_Sub2_slope = plot_dipeptides_with_trendline(
        dipeptides2, time_points2, df_Captopril2, df_GlyGly2, df_Buff_Ace_Sub2)

    concentrations_map = {'A': 2187 / 10e6, 'B': 729 / 10e6, 'C': 243 / 10e6, 'D': 81 / 10e6, 'E': 27 / 10e6,
                          'F': 9 / 10e6, 'G': 3 / 10e6, 'H': 1 / 10e6}

    slopes = [('LeuIle', df_LeuIle_slopes), ('LeuLeu', df_LeuLeu_slopes), ('LeuVal', df_LeuVal_slopes),
              ('MetVal', df_MetVal_slopes), ('MetIle', df_MetIle_slopes), ('MetLeu', df_MetLeu_slopes)]

    # plot the slopes of the linear fits against the concentrations and fit with sigmoid. use the same color for the same dipeptide
    # x axis (logarithmic) is the concentration, y axis is the slope
    # Setup a 3x2 grid of subplots
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 12), sharey=True, sharex=True)
    axes = axes.flatten()  # Flatten the array to easily iterate over it
    colors = ['red', 'blue', 'purple', 'green', 'pink', 'orange']

    for i, (name, df) in enumerate(slopes):
        if i < 3:
            captopril_slope = Captopril1_slope
            glygly_slope = GlyGly1_slope
            buffacesub_slope = Buff_Ace_Sub1_slope
        else:
            captopril_slope = Captopril2_slope
            glygly_slope = GlyGly2_slope
            buffacesub_slope = Buff_Ace_Sub2_slope

        ax = axes[i]  # Select the current Axes object
        color = colors[i]
        y_slopes = np.array(df[0])
        x_concentrations = np.array([concentrations_map[conc] for conc in df.index])

        ax.plot(x_concentrations, y_slopes, 'x', label=name, color=color)
        fit_sigmoid(ax, color, concentrations=x_concentrations, slopes=y_slopes)

        # plot a line at the control slopes
        ax.axhline(y=captopril_slope, color='red', linestyle='--', label='Captopril')
        ax.axhline(y=glygly_slope, color='blue', linestyle='--', label='GlyGly')
        ax.axhline(y=buffacesub_slope, color='darkblue', linestyle='--', label='Buff_Ace_Sub')
        # draw a line at 50% of Buff_Ace_Sub
        ax.axhline(y=buffacesub_slope / 2, color='darkblue', linestyle='--', label='Buff_Ace_Sub/2')

        # do a linear fit of four highest concentrations (A-D) to get the trendline of the linear part
        x = np.log10(x_concentrations)
        y = y_slopes
        params, _ = curve_fit(linear, x[:4], y[:4])
        IC50 = (buffacesub_slope / 2 - params[1]) / params[0]
        # draw this line until it crosses the Buff_Ace_Sub/2 line and label it with "IC50=... µM" (insert that value)
        x_fit = np.linspace(min(x), max(x), 100)

        trendline = linear(x_fit, *params)
        ax.plot(10 ** x_fit, trendline, '--', color=color)
        # calculate the IC50 value
        IC50 = (10 ** IC50)  # convert to micromolar

        ax.axvline(x=IC50, color=color, linestyle='--', label='IC50 = {:.2f} µM'.format(IC50 * 10e6))

        ax.set_xscale('log')
        ax.set_xticks([1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3])
        ax.set_xlabel('Concentration [µM]')
        ax.set_ylabel('Slope')
        ax.legend()

    # Adjust layout
    plt.tight_layout()
    plt.savefig("slopes.png")
    plt.show()
