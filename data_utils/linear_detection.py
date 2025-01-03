import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

#Detects linear behavior by fitting linear model with given r2 and slope tolerance
#Returns list of tuples: (start_index, end_index, slope, r2) for linear segments
def sliding_window_linear_detection(x, y, window_size=100, r2_threshold=0.95, slope_tolerance=0.1):
    
    linear_segments = []
    slopes = []
    
    n = len(x)
    
    if window_size >= n:
        raise ValueError("Window size must be smaller than the length of the data.")
    
    # Sliding window loop
    i = 0
    while i + window_size <= n:
        #Reshape x values for input in Linear_Regression function
        x_window = np.array(x[i:i+window_size]).reshape(-1, 1)
        y_window = np.array(y[i:i+window_size])
        
        # Fit linear regression
        model = LinearRegression().fit(x_window, y_window)
        y_pred = model.predict(x_window)
        
        # Calculate R2 value and slope
        r2 = r2_score(y_window, y_pred)
        slope = model.coef_[0]
        
        # Check for linearity 
        if r2 >= r2_threshold:
            if not slopes or abs(slope - slopes[-1]) <= slope_tolerance:
                linear_segments.append((i, i+window_size, slope, r2))
        
        slopes.append(slope)
        i += window_size // 2  
    
    return linear_segments

#Plots original data (blue) and linear segments (red)
def plot_sliding_window_linear_detection(x, y, linear_segments):
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label='Original Data', color='blue')
    
    for (start, end, slope, r2) in linear_segments:
        plt.plot(x[start:end], y[start:end], color='red', label=f'Linear Region (R²={r2:.2f})' if start == 0 else "")
    
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.title('Sliding Window Linear Detection')
    plt.show()
    plt.savefig('linear_segments.png')

#Finds optimal window size to find linear segments
def find_optimal_window_size(x, y, window_sizes, r2_threshold=0.95, slope_tolerance=0.1):
    best_window_size = None
    best_r2 = -np.inf

    for window_size in window_sizes:
        # Run sliding window linear detection
        linear_segments = sliding_window_linear_detection(x, y, window_size=window_size,
                                                          r2_threshold=r2_threshold,
                                                          slope_tolerance=slope_tolerance)
        
        if linear_segments:
            # Compute average R² value for the found segments
            average_r2 = np.mean([segment[3] for segment in linear_segments])
            
            if average_r2 > best_r2:
                best_r2 = average_r2
                best_window_size = window_size
        else:
            # No valid segments found
            average_r2 = -np.inf

    return best_window_size


def steps_to_reach_trend(x,y, slope_tolerance =0.1, r2_tolerance = 0.1):
    points_needed = 0
    melting_point_index = 0
 
    for i in range(len(x)):
        if y[i] >= 2673: #Find melting point index
            melting_point_index = i
            break
    x_global = np.array(x[0:melting_point_index]).reshape(-1, 1)
    y_global= np.array(y[0:melting_point_index])
        
    # Fit linear regression
    model_global = LinearRegression().fit(x_global, y_global)
    y_pred_global = model_global.predict(x_global)
    r2_global = r2_score(y_global, y_pred_global)
    slope_global = model_global.coef_[0]
    
    for i in range(2, melting_point_index):
        x_local = np.array(x[0:i]).reshape(-1,1)
        y_local = np.array(y[0:i])
        model_local = LinearRegression().fit(x_local, y_local)
        y_pred_local = model_local.predict(x_local)
        r2_local = r2_score(y_local, y_pred_local)
        slope_local = model_local.coef_[0]

        if abs(slope_global - slope_local) <= slope_tolerance and abs(r2_global-r2_local) <= r2_tolerance:
            points_needed = i
            break
    steps_needed = x[points_needed]-x[0]
    return steps_needed, slope_global

