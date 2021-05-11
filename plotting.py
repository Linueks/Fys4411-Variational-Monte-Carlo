import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def look_through(filename):
    with open(filename) as infile:
        lines = infile.readlines()
        times = np.zeros(len(lines))
        accept_ratios = np.zeros(len(lines))
        alphas = np.zeros(len(lines))
        energies = np.zeros(len(lines))
        gradients = np.zeros(len(lines))

        for k, line in enumerate(lines):
            variables = line.split(',')
            times[k] = variables[0].split(':')[1].rstrip('ms')
            #print(times[k])
            accept_ratios[k] = variables[5].split(':')[1]
            alphas[k] = (variables[7].split(':')[1])
            energies[k] = (variables[8].split(':')[1])
            gradients[k] = variables[9].split(':')[1]

    return times, accept_ratios, alphas, energies, gradients




def plot_comparison():
    local_energy_plot, ax = plt.subplots(figsize=(8, 5))
    sim_time_plot, ax2 = plt.subplots(figsize=(8, 5))
    gradient_plot, ax3 = plt.subplots(figsize=(8, 5))


    line_styles = ['solid', 'dashed', 'dashdot']


    sim_time_avg = np.zeros(shape=(len(dims), len(ns)))
    sim_time_std = np.zeros_like(sim_time_avg)
    confidence_interval = 1.95

    print(sim_time_avg)

    for i, d in enumerate(dims):
        for j, n in enumerate(ns):
                times, accept_ratios, alphas, energies, gradients = look_through(f'data\\{d}D_N{n}_analytic_brute')
                sim_time_avg[i, j] = np.mean(times)
                sim_time_std[i, j] = np.std(times)

                ax.plot(alphas, energies, label=f'd = {d}, n = {n}',
                        linestyle=line_styles[d - 1])

                ax3.plot(alphas, gradients, label=f'd = {d}, n = {n}',
                        linestyle=line_styles[d - 1])


        print(f'd: {d}, sim time average: ', sim_time_avg[i, :])
        print(f'd: {d}, sim time deviate: ', sim_time_std[i, :])
        ax2.plot(ns, sim_time_avg[i, :], label=f'd = {d}')
        ax2.fill_between(ns,
            sim_time_avg[i, :] - confidence_interval * sim_time_avg[i, :] / sim_time_std[i, :],
            sim_time_avg[i, :] + confidence_interval * sim_time_avg[i, :] / sim_time_std[i, :], alpha=0.1)


    ax.set_xlabel('$\\alpha$')
    ax.set_ylabel('$E_L$')
    ax.set_title('Local Energy as a Function of $\\alpha$ \n with Numeric Derivative')
    ax.legend(loc='upper right')

    ax2.set_xlabel('N')
    ax2.set_ylabel('Simulation Time [ms]')
    ax2.set_title('Simulation Time with Numeric Derivative \n as a Function of Dimension and Particle Number')
    ax2.legend(loc='upper right')

    ax3.set_xlabel('$\\alpha$')
    ax3.set_ylabel('$\\frac{\\partial E_L}{\\partial \\alpha}$')
    ax3.set_title('Energy Gradient at End of Simulation')
    ax3.legend(loc='upper right')

    plt.show()



def plot_gradient_descent():
    line_styles = ['solid', 'dashed', 'dashdot']
    dims = np.array([1, 2, 3])
    ns = np.array([1, 10, 100, 500])


    gradient_plot, ax = plt.subplots(figsize=(8, 5))
    sim_time_avg = np.zeros(shape=(len(dims), len(ns)))
    sim_time_std = np.zeros_like(sim_time_avg)



    for i, d in enumerate(dims):
        print(d)
        for j, n in enumerate(ns):
            print(n)
            times, accept_ratios, alphas, energies, gradients = look_through(f'data\\{d}D_N{n}_steepest')
            sim_time_avg[i, j] = np.mean(times)
            sim_time_std[i, j] = np.std(times)
            print(gradients)
            ax.plot(range(len(gradients)),
                    gradients,
                    label = f'd = {d}, n = {n}',
                    linestyle=line_styles[d - 1])

        print(f'd: {d}, sim time average: ', sim_time_avg[i, :])
        print(f'd: {d}, sim time deviate: ', sim_time_std[i, :])



    ax.set_xlabel('Iterations')
    ax.set_ylabel('$\\frac{\\partial E_L}{\\partial \\alpha}$')
    ax.set_title('Convergence of Local Energy Gradient')
    ax.legend()
    plt.show()






if __name__=='__main__':

    #plot_comparison()
    plot_gradient_descent()
