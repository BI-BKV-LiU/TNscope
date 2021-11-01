#!/usr/bin/awk -f

BEGIN{
    OFS="\t";
}
{
	chr = $1
	start = $2
	end = $3
	id = $4
	symbol = $5

	region[chr"_"start"_"end] = symbol

}
END{
	for (i in region) {
		split(i,position,"_")
		print position[1], position[2], position[3], region[i]
	}
}