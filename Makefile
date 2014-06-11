doc:
	doxygen

publish: doc
	cd html
	git add -A
	git commit -m "Automatic publishing"
	git push

clean:
	rm -rf html/*

reinit:
	rm -rf html
	git clone -b gh-pages git@github.com:tonigi/vmd_extensions.git html
