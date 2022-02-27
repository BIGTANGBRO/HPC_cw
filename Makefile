.PHONY: clean All

All:
	@echo "----------Building project:[ CW - Debug ]----------"
	@"$(MAKE)" -f  "CW.mk"
clean:
	@echo "----------Cleaning project:[ CW - Debug ]----------"
	@"$(MAKE)" -f  "CW.mk" clean
