import requests
import pandas as pd
import sys

class RowData:
    def __init__(self, rank, score, miRNA, gene, description):
        self.rank = rank
        self.score = score
        self.miRNA = miRNA
        self.gene = gene
        self.description = description


def getBlock(text, tag, index):
    # Return the HTML block matching the tag argument. If more than one exists,
    # we get the one at the argument index.
    currentIndex = 0
    tagWithoutBrackets = tag[1:len(tag) - 1]
    for i in range(0, len(text) - len(tag)):
        if text[i] == '<':
            textFromHereToEnd = text[i: len(text)]
            currentTagWithBrackets = textFromHereToEnd[:textFromHereToEnd.index(">") + 1]
            currentTagWithoutBrackets = currentTagWithBrackets[1:len(currentTagWithBrackets) - 1]
            if currentTagWithoutBrackets.startswith(tagWithoutBrackets):
                if currentIndex == index:
                    closingTag = '</' + tagWithoutBrackets + '>'
                    closingTagIndex = textFromHereToEnd.find(closingTag)
                    blockWithBrackets = textFromHereToEnd[:closingTagIndex + len(closingTag)]
                    return blockWithBrackets
                currentIndex += 1


def getRawDataFromBlock(htmlBlock):
    # Get the actual raw data encapsulated an HTML block. Raw data is surrounded
    # by nested pairs of opening and closing tags.
    while True:
        openingBracketIndex = htmlBlock.find('<')
        closingBracketIndex = htmlBlock.find('>')
        if openingBracketIndex == -1 and closingBracketIndex == -1:
            return htmlBlock
        htmlBlock = htmlBlock[:openingBracketIndex] + htmlBlock[closingBracketIndex + 1:]


def getData(address):
    # Send the HTTP request, receive the response, parse the output, extract
    # relevant information and enclose it in a pandas data frame.
    r = requests.get(address)
    tableBlock = getBlock(r.text, '<table>', 1)
    dataFrame = pd.DataFrame({
        "rank": [],
        "score": [],
        "miRNA": [],
        "gene": [],
        "description": []
    })
    rowIndex = 0
    while True:
        row = getBlock(tableBlock, '<tr>', rowIndex)
        if row is None:
            break
        rowIndex += 1
        if rowIndex == 1:
            continue
        rankColumn = getBlock(row, '<td>', 1)
        scoreColumn = getBlock(row, '<td>', 2)
        miRNAColumn = getBlock(row, '<td>', 3)
        geneColumn = getBlock(row, '<td>', 4)
        descriptionColumn = getBlock(row, '<td>', 5)

        rank = getRawDataFromBlock(rankColumn)
        score = getRawDataFromBlock(scoreColumn)
        miRNA = getRawDataFromBlock(miRNAColumn)
        gene = getRawDataFromBlock(geneColumn)
        description = getRawDataFromBlock(descriptionColumn)
        rowDataFrame = pd.DataFrame({
            "rank": [rank],
            "score": [score],
            "miRNA": [miRNA],
            "gene": [gene],
            "description": [description]
        })
        dataFrame = pd.concat([dataFrame, rowDataFrame])
    return dataFrame


if __name__ == '__main__':
    # Replace the following with the miRNA name you want
    #miRNA = 'rno-miR-137-5p'
    # The script receives the miRNA as a single argument
    miRNA = sys.argv[1]
    data = getData(
        'http://mirdb.org/cgi-bin/search.cgi?'
        'searchType=miRNA'
        '&searchBox=' + miRNA +
        '&full=1'
        '&fbclid=IwAR0ZX67snfY2UW_r6H6GklI3lQ5gC2Jy_KrUVW3kvriaVJfbPqSy0ph4SGA')
    #data.to_csv('tab_separated_data.csv', sep='\t')
    data.to_csv(miRNA + ".tab", sep='\t')

