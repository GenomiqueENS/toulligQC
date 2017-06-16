"""
Creates a report in the form of a docx file for a minion run
"""

import docx
import os
import parser


def docxs(selection_barcode,date, flowcell_id, barcode_present):
    dico_path = parser.file_path_initialization()
    report_writing_directory = dico_path['result_directory']
    design_file_directory = dico_path['design_file_directory']
    ##### Watch out change in the directory#####
    os.chdir(report_writing_directory)
    doc = docx.Document()
    doc.add_heading('Rapport run Minion\n', level=1)
    doc.add_paragraph('Flowcell id: '+flowcell_id)
    doc.add_paragraph('Date: '+date)
    doc.add_paragraph('Histogram of read counts according to the type read')
    doc.add_picture('images/image1.png', width=docx.shared.Inches(6),height=docx.shared.Cm(9))
    doc.add_paragraph('Phred score according to the read type')
    doc.add_picture('images/image2.png', width=docx.shared.Inches(5),height=docx.shared.Cm(9))
    doc.add_paragraph('Channel counts')
    doc.add_picture('images/image3.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))
    doc.add_paragraph('Curve representing the reads number produced against the time')
    doc.add_picture('images/image4.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))
    doc.add_paragraph('Representation of the channels occupation')
    doc.add_picture('images/image6.png', width=docx.shared.Inches(8), height=docx.shared.Cm(11))
    doc.add_paragraph('Layout')
    doc.add_picture(design_file_directory+'layout.png', width=docx.shared.Inches(6),height=docx.shared.Cm(9))
    doc.add_paragraph('Read size')
    doc.add_picture('images/image7.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))

    if barcode_present:
        doc.add_paragraph("Pie chart representing barcodes")
        doc.add_picture('images/image5.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))
        selection_barcode = selection_barcode[:-1]
        doc.add_paragraph('Read length boxplot')
        for barcode_image in selection_barcode:
            doc.add_picture('images/image_{}.png'.format(barcode_image), width=docx.shared.Inches(5), height=docx.shared.Cm(9))

# height=docx.shared.Cm(4))
        doc.paragraphs[0].runs[0].underline = True
        doc.paragraphs[3].runs[0].underline = True
        doc.paragraphs[5].runs[0].underline = True
        doc.paragraphs[7].runs[0].underline = True
        doc.paragraphs[9].runs[0].underline = True
        doc.paragraphs[11].runs[0].underline = True
        doc.paragraphs[13].runs[0].underline = True
        doc.paragraphs[15].runs[0].underline = True
        doc.paragraphs[17].runs[0].underline = True
    doc.save('Rapport.docx')

