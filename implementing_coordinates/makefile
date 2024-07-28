THEME=serif

all: presentation.html

presentation.html: presentation.md
	pandoc \
	--slide-level 2 \
	-t revealjs \
	-s \
	-o presentation.html \
	presentation.md \
	-V revealjs-url=https://unpkg.com/reveal.js \
	--include-in-header=slides.css \
	-V theme=$(THEME)

clean:
	rm -f *.html