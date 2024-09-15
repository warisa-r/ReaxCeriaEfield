import numpy as np
import matplotlib.pyplot as plt

def single_plot(x_array, y_array, x_label, y_label, title, file_name):
    plt.figure(figsize=(8, 6))
    plt.plot(x_array, y_array)

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    plt.grid(True)

    #plt.savefig(file_name, format=file_name.split('.')[-1], dpi=300)
    plt.show()



def double_plot(x_array1, y_array1, x_label1, y_label1, title1,
                      x_array2, y_array2, x_label2, y_label2, title2, 
                      file_name):
    # Create subplots with 1 row and 2 columns
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    axs[0].plot(x_array1, y_array1)
    axs[0].set_xlabel(x_label1)
    axs[0].set_ylabel(y_label1)
    axs[0].set_title(title1)
    axs[0].grid(True)

    axs[1].plot(x_array2, y_array2)
    axs[1].set_xlabel(x_label2)
    axs[1].set_ylabel(y_label2)
    axs[1].set_title(title2)
    axs[1].grid(True)
    plt.tight_layout()

    #plt.savefig(file_name, format=file_name.split('.')[-1], dpi=300)

    plt.show()


def triple_plot(x_array1, y_array1, x_label1, y_label1, title1,
                x_array2, y_array2, x_label2, y_label2, title2,
                x_array3, y_array3, x_label3, y_label3, title3,
                file_name):
    # Create subplots with 1 row and 3 columns
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    # Plot on the first subplot
    axs[0].plot(x_array1, y_array1)
    axs[0].set_xlabel(x_label1)
    axs[0].set_ylabel(y_label1)
    axs[0].set_title(title1)
    axs[0].grid(True)

    # Plot on the second subplot
    axs[1].plot(x_array2, y_array2)
    axs[1].set_xlabel(x_label2)
    axs[1].set_ylabel(y_label2)
    axs[1].set_title(title2)
    axs[1].grid(True)

    # Plot on the third subplot
    axs[2].plot(x_array3, y_array3)
    axs[2].set_xlabel(x_label3)
    axs[2].set_ylabel(y_label3)
    axs[2].set_title(title3)
    axs[2].grid(True)

    plt.tight_layout()

    #plt.savefig(file_name, format=file_name.split('.')[-1], dpi=300)

    plt.show()