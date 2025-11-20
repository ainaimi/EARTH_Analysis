# Makefile for EARTH Analysis Pipeline
# Run entire pipeline with: make all
# Clean outputs with: make clean

# Variables
R := Rscript
QUARTO := quarto render code/EARTH_Analysis_Report.qmd --output-dir .

# Directories
CODE_DIR := code
DATA_DIR := data
FIG_DIR := figures
MISC_DIR := misc

# Data outputs from each step
DATA_OUTPUTS := $(DATA_DIR)/afc_clean_trunc.Rdata $(DATA_DIR)/afc_clean_notrunc.Rdata
CHEM_OUTPUTS := $(FIG_DIR)/.chem_analysis_done
IF_OUTPUTS := $(MISC_DIR)/fit_mu.RDS $(MISC_DIR)/fit_pi.RDS
ANALYSIS_OUTPUTS := $(MISC_DIR)/.analysis_done

# Final report
REPORT := EARTH_Analysis_Report.html

# Default target: run entire pipeline and generate report
.PHONY: all
all: $(REPORT)
	@echo "======================================"
	@echo "Pipeline completed successfully!"
	@echo "Report available at: $(REPORT)"
	@echo "======================================"

# Step 1: Data management
$(DATA_OUTPUTS): $(CODE_DIR)/1.\ AFC_data_man.R
	@echo "======================================"
	@echo "Step 1: Data Management"
	@echo "======================================"
	cd $(CODE_DIR) && $(R) "1. AFC_data_man.R"
	@touch $(DATA_OUTPUTS)

# Step 2: Chemical analysis
$(CHEM_OUTPUTS): $(CODE_DIR)/2.\ AFC_chem_analysis.R $(DATA_OUTPUTS)
	@echo "======================================"
	@echo "Step 2: Chemical Analysis"
	@echo "======================================"
	cd $(CODE_DIR) && $(R) "2. AFC_chem_analysis.R"
	@mkdir -p $(FIG_DIR)
	@touch $(CHEM_OUTPUTS)

# Step 3: IF score generation (this may take a while)
$(IF_OUTPUTS): $(CODE_DIR)/3.\ AFC_IF_scores_gen.R $(DATA_OUTPUTS)
	@echo "======================================"
	@echo "Step 3: IF Score Generation"
	@echo "This may take several minutes..."
	@echo "======================================"
	cd $(CODE_DIR) && $(R) "3. AFC_IF_scores_gen.R"
	@touch $(IF_OUTPUTS)

# Step 4: Heterogeneity analysis
$(ANALYSIS_OUTPUTS): $(CODE_DIR)/4.\ AFC_IF_scores_analysis.R $(IF_OUTPUTS)
	@echo "======================================"
	@echo "Step 4: Heterogeneity Analysis"
	@echo "======================================"
	cd $(CODE_DIR) && $(R) "4. AFC_IF_scores_analysis.R"
	@mkdir -p $(MISC_DIR)
	@touch $(ANALYSIS_OUTPUTS)

# Generate final report
$(REPORT): $(CODE_DIR)/EARTH_Analysis_Report.qmd $(ANALYSIS_OUTPUTS)
	@echo "======================================"
	@echo "Generating Final Report"
	@echo "======================================"
	$(QUARTO)
	@echo "Report generated: $(REPORT)"

# Individual step targets (can run steps independently)
.PHONY: step1 step2 step3 step4 report

step1: $(DATA_OUTPUTS)

step2: $(CHEM_OUTPUTS)

step3: $(IF_OUTPUTS)

step4: $(ANALYSIS_OUTPUTS)

report: $(REPORT)

# Clean outputs (be careful with this!)
.PHONY: clean clean-all clean-report

clean-report:
	@echo "Removing report..."
	rm -f $(REPORT)
	rm -f $(CODE_DIR)/EARTH_Analysis_Report.html
	rm -rf $(CODE_DIR)/EARTH_Analysis_Report_files

clean:
	@echo "Removing intermediate marker files..."
	rm -f $(FIG_DIR)/.chem_analysis_done
	rm -f $(MISC_DIR)/.analysis_done

clean-all: clean clean-report
	@echo "WARNING: This will remove ALL generated outputs!"
	@echo "Press Ctrl+C to cancel, or wait 3 seconds to proceed..."
	@sleep 3
	@echo "Removing all generated files..."
	rm -rf $(FIG_DIR)/*
	rm -f $(MISC_DIR)/*.RDS
	rm -f $(DATA_DIR)/*.Rdata
	@echo "Clean complete. Data, figures, and outputs removed."

# Help target
.PHONY: help
help:
	@echo "EARTH Analysis Pipeline - Makefile Commands"
	@echo "==========================================="
	@echo ""
	@echo "Main targets:"
	@echo "  make all          - Run entire pipeline and generate report"
	@echo "  make report       - Generate HTML report (assumes steps 1-4 complete)"
	@echo ""
	@echo "Individual steps:"
	@echo "  make step1        - Data management"
	@echo "  make step2        - Chemical analysis and outlier detection"
	@echo "  make step3        - IF score generation (slow)"
	@echo "  make step4        - Heterogeneity analysis and plots"
	@echo ""
	@echo "Cleaning:"
	@echo "  make clean        - Remove intermediate marker files"
	@echo "  make clean-report - Remove report only"
	@echo "  make clean-all    - Remove ALL generated outputs (use with caution!)"
	@echo ""
	@echo "Other:"
	@echo "  make help         - Show this help message"
