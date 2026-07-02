import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plot_observability_windows():
    # Load the CSV file, parsing the UTC strings into datetime objects
    df = pd.read_csv('voyager_station_observability_windows.csv', parse_dates=['start_utc', 'end_utc'])
    
    fig, ax = plt.subplots(figsize=(12, 5))
    fig.canvas.manager.set_window_title('DSN Tracking Passes')
    
    # Get unique stations and sort them for a consistent Y-axis
    stations = sorted(df['station'].unique(), reverse=True)
    colors = plt.cm.tab10.colors

    # Plot each tracking pass as a horizontal line
    for idx, row in df.iterrows():
        y_idx = stations.index(row['station'])
        ax.hlines(y=y_idx, xmin=row['start_utc'], xmax=row['end_utc'], 
                  linewidth=25, color=colors[y_idx % len(colors)], capstyle='round')

    # Formatting the Y-axis
    ax.set_yticks(range(len(stations)))
    ax.set_yticklabels(stations, fontsize=12, fontweight='bold')
    
    # Formatting the X-axis
    ax.set_xlabel('UTC Time', fontsize=11)
    ax.set_title('Voyager 1 DSN Station Observability Windows', fontsize=14)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M'))
    fig.autofmt_xdate() # Auto-rotates the dates so they don't overlap
    
    # Add a grid for easier visual alignment of overlaps
    ax.grid(True, axis='x', linestyle='--', alpha=0.7)
    
    # Add vertical padding
    ax.set_ylim(-0.5, len(stations) - 0.5)
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot_observability_windows()