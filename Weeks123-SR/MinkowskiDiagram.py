import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from InvariantInterval import interval
from LorentzBoost import lorentz_boost, lorentz_factor

mpl.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "axes.labelsize": 14,
    "font.size": 14,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

import random
 
def random_colour_generator():
    r = 0
    g = random.randint(0, 125)/255
    b = random.randint(0, 255)/255
    return (r, g, b)

def light_cone_plot(events, size, shade_common_future = False, shade_common_past = False):
    fig, ax = plt.subplots(figsize = [6,6])

    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)

    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    ax.set_xlabel('x', loc='right', fontname = 'DejaVu Sans')
    ax.set_ylabel('t', rotation=0, loc='top', fontname = 'DejaVu Sans')

    c = np.linspace(-size, size, 100)
    # origin light cone
    # ax.plot(c, c, linestyle = '--', alpha = 0.5, color = 'blue')
    # ax.plot(c, -c, linestyle = '--', alpha = 0.5, color = 'blue')

    for event in events:
        x_coord = event[1]
        t_coord = event[0]

        colour = random_colour_generator()
        ax.scatter(x_coord, t_coord, marker = 'o', s = 10, color = colour, zorder = 11, label = 'Event ' + str(events.index(event) + 1))
        ax.plot(c, t_coord + c - x_coord, linestyle = '--', color = colour, zorder = 10, label = 'Light cone of event ' + str(events.index(event) + 1))
        ax.plot(c, t_coord - c + x_coord, linestyle = '--', color = colour, zorder = 10)

    # if shade_common_future and len(events) > 1:
    #     t_min = max(t for t, _ in events)
    #     t_vals = np.linspace(t_min, size, 500)

    #     left = np.full_like(t_vals, -np.inf)
    #     right = np.full_like(t_vals, np.inf)

    #     for t_i, x_i in events:
    #     # Only defined for t_vals ≥ t_i
    #         mask = t_vals >= t_i
    #         left_event = np.full_like(t_vals, -np.inf)
    #         right_event = np.full_like(t_vals, np.inf)

    #         left_event[mask] = x_i - (t_vals[mask] - t_i)
    #         right_event[mask] = x_i + (t_vals[mask] - t_i)

    #         left = np.maximum(left, left_event)
    #         right = np.minimum(right, right_event)

    #     valid = left < right
    #     ax.fill_between(t_vals[valid], left[valid], right[valid], color='orange', alpha=0.3)


    # if shade_common_past and len(events) > 1:
    #     t_max = min(t for t, _ in events)
    #     t_vals = np.linspace(-size, t_max, 500)

    #     left = np.full_like(t_vals, -np.inf)
    #     right = np.full_like(t_vals, np.inf)

    #     for t_i, x_i in events:
    #     # Only defined for t_vals ≤ t_i
    #         mask = t_vals <= t_i
    #         left_event = np.full_like(t_vals, -np.inf)
    #         right_event = np.full_like(t_vals, np.inf)

    #         left_event[mask] = x_i - (t_i - t_vals[mask])
    #         right_event[mask] = x_i + (t_i - t_vals[mask])

    #         left = np.maximum(left, left_event)
    #         right = np.minimum(right, right_event)

    #     valid = left < right
    #     ax.fill_between(t_vals[valid], left[valid], right[valid], color='purple', alpha=0.3)

    grid_ticks = range(-size, size + 1, 1)
    ax.set_xticks(grid_ticks, minor=False)
    ax.set_yticks(grid_ticks, minor=False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    ax.grid(True, which='both', axis='both', linestyle='--', alpha=0.3)
    ax.text(-0.1, -0.1, '0', fontsize=12, ha='right', va='top', fontname='DejaVu Sans')
    ax.legend(loc = 'best', fontsize = 8)
    ax.figure.canvas.draw_idle()

    plt.show()

test = [[0, 1], [1, 4]]
# light_cone_plot(test, 10, True, True)


def minkowski_diagram(events, vs, frame, size, hyperbolae = False, event_coordinates = False):
    fig, ax = plt.subplots(figsize = [6,6])


    step = 2  # or 1 or any appropriate value
    xticks = np.arange(-size, size + 1, step)
    yticks = np.arange(-size, size + 1, step)
    xtick_labels = [str(x) if x != 0 else '' for x in xticks]
    ytick_labels = [str(y) if y != 0 else '' for y in yticks]
    plt.xticks(xticks, xtick_labels)
    plt.yticks(yticks, ytick_labels)
    ax.set_xlim(-size, size)
    ax.set_ylim(-size, size)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_xlabel(r'$x_0$', loc='right', fontname = 'DejaVu Sans')
    ax.set_ylabel(r'$t_0$', rotation=0, loc='top', fontname = 'DejaVu Sans')
    plt.text(-0.4, -0.4, 'O', ha='center', va='center', fontsize=12, zorder = 15, weight = 'bold')

    c = np.linspace(-size, size, 100)

    for n in range(len(events)):
        event = events[n]
        x_coord = event[1]
        t_coord = event[0]
        colour = random_colour_generator()
        ax.scatter(x_coord, t_coord, marker = 'o', s = 15, color = colour, zorder = 11, label = 'Event ' + str(events.index(event) + 1))
        # ax.plot(c, t_coord + c - x_coord, linestyle = '--', color = colour, zorder = 10, label = 'Light cone of event ' + str(events.index(event) + 1))
        # ax.plot(c, t_coord - c + x_coord, linestyle = '--', color = colour, zorder = 10)

    rest_coords = []
    for event in events:
        rest_coords.append(event)
    
    all_coords = [rest_coords]
    
    for l in range(len(vs)):
        v = vs[l]
        if np.abs(v) <= 1:
            plt.plot(c, c/v, linestyle = ':', color = 'k', zorder = 15) # t' axis has x' = 0
            plt.plot(c, v*c, linestyle = ':', color = 'k', zorder = 15) # x' axis has t' = 0
            plt.text(size, v*size, fr"$x'_{{{l + 1}}}$")
            plt.text(v*size, size, fr"$t'_{{{l + 1}}}$")

            gamma = lorentz_factor(v)
            event_coords = []
            for e in range(len(events)):
                tp = gamma*events[e][0] - gamma*v*events[e][1]
                xp = -gamma*v*events[e][0] + gamma*events[e][1]
                primed = [tp, xp]
                event_coords.append(primed)
            all_coords.append(event_coords)
        else:
            raise ValueError(f"velocity v={v} must satisfy |v| <= 1.")
    # print(all_coords)
    
    if event_coordinates == True:
        for n in range(len(events)):
            plt.text(all_coords[0][n][1], all_coords[0][n][0] - 0.2, f'({round(all_coords[frame][n][0], 2)}, {round(all_coords[frame][n][1], 2)})', ha = 'center', va = 'top', fontsize = 10, zorder = 11)
        
    if hyperbolae == True:
        plt.plot(c, np.sqrt(c**2-1), color = 'k', linestyle = '--', label = r'Calibration hyperbolae for an interval of $\pm 1$')
        plt.plot(c, -np.sqrt(c**2-1), color = 'k', linestyle = '--')
        plt.plot(c, np.sqrt(c**2+1), color = 'k', linestyle = '--')
        plt.plot(c, -np.sqrt(c**2+1), color = 'k', linestyle = '--')
        # plt.scatter(0, 1, s = 5, zorder = 10, color = 'k')
        # plt.scatter(1, 0, s = 5, zorder = 10, color = 'k')
        # plt.scatter(-1, 0, s = 5, zorder = 10, color = 'k')
        # plt.scatter(0, -1, s = 5, zorder = 10, color = 'k')

    plt.legend(loc = 'best', fontsize = 8 )
    plt.show()

    print(all_coords)
    print(len(all_coords))

tester = [[0, 4], [3, 1], [-3, -2]]

tester_v = [0.5, 0.3]
minkowski_diagram(tester, tester_v, 2, 8, False, True)
        
#AXIS LINES CONNECTING TO POINTS
        




