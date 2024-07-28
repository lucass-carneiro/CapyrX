THEME=dracula

all: presentation.html

presentation.html: presentation.md slides.css
	pandoc \
	--slide-level 2 \
	-t revealjs \
	-s \
	-o presentation.html \
	presentation.md \
	--include-in-header=slides.css

clean:
	rm -f *.html
