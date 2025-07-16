doc:
	doxygen

# Needs doxygen 1.8.13. TCL support was removed in 1.8.18. 1.8.17 should be the last.

# 1.8.17 not in conda
#    .16 empty page
#    .15 not in conda
#    .14 works
#    .13 works

clean:
	rm -rf html/*
