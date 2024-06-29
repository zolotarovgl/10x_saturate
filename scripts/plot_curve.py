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
    
    return pred_train, pred_new, popt

# Argument parsing
parser = argparse.ArgumentParser(description='Fit a saturation model and plot the saturation curve.')
parser.add_argument('infile', help='Input file containing data (CSV format)')
parser.add_argument('plotfile', help='Output plot file (PNG format)')
parser.add_argument('--maxx', type=float, default=5, help='Factor by which to increase coverage')
parser.add_argument('--n_points', type=int, default=200, help='Number of points to predict saturation')
parser.add_argument('--target', type=float, default=0.7, help='Target saturation')
args = parser.parse_args()

# Read input data
try:
    d = pd.read_csv(args.infile, sep='\t')
except FileNotFoundError:
    print(f"Error: File '{args.infile}' not found.")
    exit(1)

# Predict using fit_model function
da = pd.DataFrame({'ninput': np.linspace(min(d['ninput']), max(d['ninput']) * args.maxx, num=args.n_points)})
pred_train, pred_new, popt = fit_model(d, da)

highest_saturation = np.max(pred_new)

# Calculate R-squared and RMSE for training data
r2_train = r2_score(d['sat'], pred_train)
rmse_train = np.sqrt(mean_squared_error(d['sat'], pred_train))

# Calculate R-squared for new data (testing data)
r2_new = r2_score(d['sat'], pred_train)  # Same as training data
rmse_new = np.sqrt(mean_squared_error(d['sat'], pred_train))  # Same as training data


# Plotting with metrics
plt.figure(figsize=(10, 6))
plt.plot(da['ninput'] / 1e6, pred_new, color='gray', label='Predicted Saturation (New Data)', alpha=0.5)
plt.scatter(d['ninput'] / 1e6, pred_train, color='red', label='Predicted Saturation (Training Data)')
plt.xlabel('Coverage, M reads')
plt.ylabel('Saturation')
plt.axvline(x=d.iloc[-1]['ninput'] / 1e6, color='red', linestyle='--', label='Last Data Point')
plt.axhline(y=highest_saturation, color='grey',linestyle=':', label=f'Highest Saturation: {highest_saturation:.2f}')
plt.title('Saturation Curve with Model Fit')
plt.text(0.05, 0.9, f'R-squared: {r2_train:.3f}\nRMSE: {rmse_train:.3f}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='top')
plt.legend()

# Save the plot
plt.savefig(args.plotfile)
plt.close()

print(f"Plot saved to {args.plotfile}")

# Report the input ninput where saturation is 0.7
def find_ninput_for_saturation(Vmax, Km, target_saturation=0.7):
    # Function to solve for ninput when saturation is target_saturation
    def saturation_eq(ninput):
        return saturation_model(ninput, Vmax, Km) - target_saturation
    
    # Initial guess for ninput (can be adjusted based on your data)
    initial_guess = 1.0
    
    # Use root finding method to find ninput
    from scipy.optimize import fsolve
    ninput_solution = fsolve(saturation_eq, initial_guess)[0]
    
    return ninput_solution


ninput_saturation = find_ninput_for_saturation(popt[0], popt[1], args.target) / 1e6
print(f"To achieve a saturation of {args.target:.2f}, ninput should be approximately: {ninput_saturation:.1f} M reads")

