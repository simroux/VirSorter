PRG = Sliding_windows_3

sliding: clean
	gcc -c $(PRG).c
	gcc -Wall $(PRG).o -o $(PRG)
	chmod +x $(PRG)

clean:
	find . -name \.o -exec rm {} \;
