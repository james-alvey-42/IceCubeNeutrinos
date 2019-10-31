import subprocess

if __name__ == '__main__':
    parameter_list = [['0.1', '0.01'], ['0.1', '0.03'], ['0.1', '0.04'], ['0.5', '0.01'], ['0.5', '0.03'], ['0.5', '0.04'], ['1.0', '0.01'], ['1.0', '0.03'], ['1.0', '0.04'], ['1.5', '0.03'], ['1.5', '0.04']]

    for params in parameter_list:
        command = 'time python include_masses_plots.py {} {}'.format(params[0], params[1])
        print(command)
        process = subprocess.Popen([command], shell=True)
