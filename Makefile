doc:
	doxygen

# Needs doxygen 1.8.13. TCL support was removed in 1.8.18. 1.8.17 should be the last.

publish: doc
	cd html; \
	git add -A; git commit -m "Automatic publishing"; git push --force

clean:
	rm -rf html/*

reinit:
	rm -rf html
	git clone -b gh-pages git@github.com:tonigi/vmd_extensions.git html
