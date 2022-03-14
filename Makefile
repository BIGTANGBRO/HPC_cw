.PHONY: clean All

All:
	@echo "----------Building project:[ coursework - Debug ]----------"
	@cd "coursework" && "$(MAKE)" -f  "coursework.mk"
clean:
	@echo "----------Cleaning project:[ coursework - Debug ]----------"
	@cd "coursework" && "$(MAKE)" -f  "coursework.mk" clean
