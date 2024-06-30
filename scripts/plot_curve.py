import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_squared_error

# Define the saturation model function
def saturation_model(ninput, Vmax, Km):
    return Vmax * ninput / (Km + ninput)

# Function to fit model and predict saturation
def fit_model(d, newdata):
    # Fit the model using curve_fit from scipy
    popt, pcov = curve_fit(saturation_model, d['ninput'], d['sat'], maxfev=10000)

    # Predict saturation for data used in training (d)
    pred_train = saturation_model(d['ninput'], *popt)

    # Predict saturation for new data (newdata)
    pred_new = saturation_model(newdata['ninput'], *popt)

    return pred_train, pred_new, popt, pcov

# Function to compute confidence intervals
def compute_confidence_intervals(x, popt, pcov, alpha=0.95):
    from scipy.stats import t
    Vmax, Km = popt
    n = len(x)
    p = len(popt)
    dof = max(0, n - p)
    tval = t.ppf((1 + alpha) / 2, dof)
    sigma = np.sqrt(np.diag(pcov))

    confs = []
    for xi in x:
        J = np.array([(xi / (Km + xi)), (-Vmax * xi / (Km + xi)**2)])
        conf = tval * np.sqrt(J @ pcov @ J.T)
        confs.append(conf)
    
    return np.array(confs)

# Argument parsing
parser = argparse.ArgumentParser(description='Fit a saturation model and plot the saturation curve.')
parser.add_argument('infile', help='Input file containing data (CSV format)')
parser.add_argument('plotfile', help='Output plot file (PNG format)')
parser.add_argument('--maxx', type=float, default=5, help='Factor by which to increase coverage')
parser.add_argument('--n_points', type=int, default=200, help='Number of points to predict saturation')
parser.add_argument('--target', type=float, default=0.7, help='Target saturation')
args = parser.parse_args()

target_saturation = args.target

# Read input data
try:
    d = pd.read_csv(args.infile, sep='\t')
except FileNotFoundError:
    print(f"Error: File '{args.infile}' not found.")
    exit(1)

# Predict using fit_model function
da = pd.DataFrame({'ninput': np.linspace(min(d['ninput']), max(d['ninput']) * args.maxx, num=args.n_points)})
pred_train, pred_new, popt, pcov = fit_model(d, da)
highest_saturation = np.max(pred_new)

# Calculate residuals
residuals = d['sat'] - pred_train

# Compute confidence intervals for the predictions
confidence_intervals = compute_confidence_intervals(da['ninput'], popt, pcov)

# Calculate R-squared and RMSE for training data
r2_train = r2_score(d['sat'], pred_train)
rmse_train = np.sqrt(mean_squared_error(d['sat'], pred_train))

# Report the input ninput where saturation is 0.7
def find_ninput_for_saturation(Vmax, Km, target_saturation=0.7):
    # Function to solve for ninput when saturation is target_saturation
    def saturation_eq(ninput):
        return saturation_model(ninput, Vmax, Km) - target_saturation

    # Initial guess for ninput (can be adjusted based on your data)
    initial_guess = 1e7  # at least 10M reads

    # Use root finding method to find ninput
    from scipy.optimize import fsolve
    ninput_solution = fsolve(saturation_eq, initial_guess)[0]

    print(f"Vmax: {Vmax:.3f}; Km: {Km:.3f}")
    print(f"Initial guess: {initial_guess}")
    print(f"ninput_solution: {ninput_solution}")

    return ninput_solution

ninput_saturation = find_ninput_for_saturation(popt[0], popt[1], target_saturation)

print(f"To achieve a saturation of {target_saturation:.2f}, ninput should be approximately: {ninput_saturation / 1e6 :.1f} M reads")

# Plotting with metrics
plt.figure(figsize=(10, 6))
plt.plot(da['ninput'] / 1e6, pred_new, color='gray', label='Projected saturation', alpha=0.5)
plt.fill_between(da['ninput'] / 1e6, pred_new - confidence_intervals, pred_new + confidence_intervals, color='gray', alpha=0.2, label='Confidence interval')
plt.scatter(d['ninput'] / 1e6, pred_train, color='red', label='Downsampled saturation')
plt.xlabel('Coverage, M reads')
plt.ylabel('Saturation')
plt.axvline(x=d.iloc[-1]['ninput'] / 1e6, color='red', linestyle=':')
plt.axvline(x=ninput_saturation / 1e6, color='blue', linestyle=':', label=f'Needed input: {ninput_saturation / 1e6:.1f}')
plt.axhline(y=target_saturation, color='lightblue', linestyle=':', label=f'Target: {target_saturation:.2f}')
plt.axhline(y=highest_saturation, color='grey', linestyle=':', label=f'Highest Saturation: {highest_saturation:.2f}')
plt.title(args.infile)
plt.text(0.05, 0.9, f'R-squared: {r2_train:.3f}\nRMSE: {rmse_train:.3f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.legend()

# Save the plot
plt.savefig(args.plotfile)
plt.close()

print(f"Plot saved to {args.plotfile}")

# Plot residuals
plt.figure(figsize=(10, 6))
plt.scatter(d['ninput'] / 1e6, residuals, color='red')
plt.axhline(0, color='gray', linestyle='--')
plt.xlabel('Coverage, M reads')
plt.ylabel('Residuals')
plt.title('Residuals Plot')
plt.savefig(args.plotfile.replace('.png', '_residuals.png'))
plt.close()

print(f"Residuals plot saved to {args.plotfile.replace('.png', '_residuals.png')}")

