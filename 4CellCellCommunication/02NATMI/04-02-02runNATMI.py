python ExtractEdges.py --species mouse --emFile /Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/ND_expressionMatrix.csv --annFile /Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/ND_metadata.csv --interDB lrc2p --coreNum 4

python ExtractEdges.py --species mouse --emFile /Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/HFD_expressionMatrix.csv --annFile /Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/HFD_metadata.csv --interDB lrc2p --coreNum 4

python DiffEdges.py \
  --refFolder /Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/ND_expressionMatrix \
  --targetFolder /Users/coellearth/Desktop/Mammary_Gland_Diet_Project/4CellCellCommunication/02NATMI/HFD_expressionMatrix \
  --interDB lrc2p \
  --weightType mean \
  --out Diff_HFD_vs_ND