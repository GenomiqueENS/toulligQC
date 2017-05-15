import docx
#'
###############directory /home/ferrato/ownCloud/Documents/fast5_1D/images###################

def docxs(selection_barcode,date, flowcell_id):
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
    doc.add_paragraph("Pie chart representing barcodes")
    doc.add_picture('images/image5.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))
    doc.add_paragraph('Curve representing the reads number produced against the time')
    doc.add_picture('images/image4.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))
    doc.add_paragraph('Representation of the channels occupation')
    doc.add_picture('images/image6.png', width=docx.shared.Inches(8), height=docx.shared.Cm(11))
    doc.add_paragraph('Layout')
    doc.add_picture('layout.png', width=docx.shared.Inches(6),height=docx.shared.Cm(9))
    doc.add_paragraph('Read size')
    doc.add_picture('images/image7.png', width=docx.shared.Inches(5), height=docx.shared.Cm(9))
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


