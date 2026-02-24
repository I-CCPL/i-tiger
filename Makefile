include settings.mk

.PHONY: all clean depend print-vars

all: $(TARGET)

$(TARGET): $(DEP_FILE) $(OBJ)
	@mkdir -p $(dir $@)
	$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

$(BUILD_DIR)/%.o: %.f90 | $(BUILD_DIR)
	$(MPIF90) $(FFLAGS) -c $< -o $@

$(BUILD_DIR):
	@mkdir -p $@

depend: $(DEP_FILE)

# Generate dependency file automatically when makedepf90 is available.
$(DEP_FILE): $(SRC) | $(BUILD_DIR)
	@set -e; \
	if command -v makedepf90 >/dev/null 2>&1; then \
	  echo "Generating $@ with makedepf90."; \
	  makedepf90 $(SRC) -nosrc -b='$(BUILD_DIR)/' > $@; \
	else \
	  echo "# makedepf90 not found; dependency generation skipped"; \
	fi

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

print-vars:
	@echo "MPIF90=$(MPIF90)"
	@echo "BUILD=$(BUILD)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "SRC=$(SRC)"
	@echo "OBJ=$(OBJ)"
	@echo "TARGET=$(TARGET)"
	@echo "DEP_FILE=$(DEP_FILE)"

-include $(DEP_FILE)
