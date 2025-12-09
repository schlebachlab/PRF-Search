for folder in "AspN_CAD" "AspN_ETD" "GluC_CAD" "GluC_ETD" "Chymo_CAD" "Chymo_ETD" "LysN_CAD" "LysN_ETD" "LysC_CAD" "LysC_ETD" "Trypsin_CAD";do
	mkdir -p ./results/"$folder"	
	rsync -av --include '*/' --include '*.tsv"' --exclude '*.pin' --exclude '*.pepXML' --exclude '*.bin' --exclude '*.png' --exclude '*.pep.xml' ./"$folder"/${folder}_2pass_results/ ./results/"$folder"/
done