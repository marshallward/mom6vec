URL=http://lab.hakim.se/reveal-js
REPO=https://github.com/hakimel/reveal.js/archive/master.zip
THEME=gfdl
FLAGS=-s \
	  -f rst -t revealjs \
	  --slide-level=2 \
	  -V revealjs-url=./reveal.js \
	  -V theme=${THEME} \
	  -V slideNumber=true \
	  --template=gfdl.revealjs \
	  --no-highlight \
	  --mathjax

#DOTFIGURES=img/gitrepos.svg img/mom_submit.svg img/mom_verify.svg
DOTFILES=$(wildcard dot/*.dot)
DOTFIGURES=$(patsubst %.dot,%.svg,$(subst dot/,img/,$(DOTFILES)))

SOURCE=$(wildcard src/*.F90) perf.stat

all: index.html reveal.js

reveal.js:
	wget -N ${REPO}
	unzip master.zip
	mv reveal.js-master reveal.js

reveal.js/css/theme/gfdl.css: gfdl.css
	cp gfdl.css reveal.js/css/theme/

index.html: slides.txt gfdl.revealjs reveal.js/css/theme/gfdl.css $(DOTFIGURES) $(SOURCE)
	pandoc ${FLAGS} $< -o $@
	sed -i 's/^" data-start-line=/"><code data-start-line=/g' index.html
	sed -i 's/^"><code>/">/g' index.html

img/%.svg: dot/%.dot
	dot -Tsvg $^ > $@

img/fixedprec.svg: dot/fixedprec.dot
	dot -Tsvg:cairo $^ > $@

clean:
	rm -f index.html 
	rm -f $(DOTFIGURES)
