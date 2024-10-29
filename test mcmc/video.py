# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 12:49:54 2024

@author: shosseini
"""




import os
import fitz  # PyMuPDF
import imageio
from PIL import Image
import numpy as np
import io  # For handling byte data
import re  # For extracting numerical values from filenames

# Define the path where your PDFs are stored
pdf_folder = r'C:\Users\shosseini\Desktop\PAPERS\RISK PAPER\test mcmc'
# Define the output file name
output_file = 'movie.mp4'

# Helper function to extract numbers from the file name for correct sorting
def extract_number(file_name):
    match = re.search(r'\d+', file_name)  # Extracts the first occurrence of a number
    return int(match.group()) if match else 0

# Get the list of PDF files in the folder, sorted numerically
pdfs = [pdf for pdf in os.listdir(pdf_folder) if pdf.endswith(".pdf")]
pdfs.sort(key=extract_number)  # Sort the files based on the numerical value in the name

# Create a list to store the image data
image_list = []

# Read each PDF and convert it to an image
for pdf_file in pdfs:
    pdf_path = os.path.join(pdf_folder, pdf_file)
    
    # Open the PDF
    pdf_document = fitz.open(pdf_path)
    
    # Iterate through pages and convert each to an image
    for page_num in range(pdf_document.page_count):
        page = pdf_document.load_page(page_num)  # Load page by page
        pix = page.get_pixmap()  # Convert the page to a pixmap (image)
        
        # Convert pixmap to byte data
        image_data = pix.tobytes(output="png")
        
        # Use io.BytesIO to handle the image byte data
        image = Image.open(io.BytesIO(image_data))  # Load image using Pillow (PIL)
        
        # Convert image to numpy array (imageio accepts this format)
        image_array = np.array(image)
        
        # Append the image data to the list
        image_list.append(image_array)

# Create the video
imageio.mimsave(output_file, image_list, fps=2)  # Adjust fps (frames per second) as needed
