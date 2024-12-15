import os
import subprocess
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit.components.v1 as components

# Sidebar for navigation
menu = st.sidebar.radio("Navigation", ["Home", "Upload Data", "Results", "Analysis", "Neo Dashboard", "Reports"])


if menu == "Home":
    st.header("Welcome to the Transcriptomics Dashboard")
    st.write("""
        This dashboard allows you to analyze and visualize transcriptomics data, including RNAseq and scRNAseq data.
        Features include:
        - **Upload your data**
        - **Perform analysis**
        - **Generate reports**
    """)
    # Placeholder for an image
    st.image("https://www.embl.org/news/wp-content/uploads/2020/05/2004_TAP-seq_Steinmetz_wordpress.jpg", use_container_width=True)

elif menu == "Upload Data":
    st.header("Upload Your Transcriptomics Data")
    file_type = st.radio("Select File Type", ["CSV/Excel", "Expression File (.h5)"])

    if file_type == "CSV/Excel":
        uploaded_file = st.file_uploader("Upload CSV or Excel file", type=["csv", "xlsx"])
        if uploaded_file:
            try:
                data = pd.read_csv(uploaded_file) if uploaded_file.name.endswith('.csv') else pd.read_excel(uploaded_file)
                st.success("File uploaded successfully!")
                st.dataframe(data.head())
            except Exception as e:
                st.error(f"Error reading the file: {e}")

    elif file_type == "Expression File (.h5)":
        uploaded_file = st.file_uploader("Upload your .h5 file", type=["h5"])
        if uploaded_file:
            file_path = f"uploaded_files/{uploaded_file.name}"
            os.makedirs("uploaded_files", exist_ok=True)

            # Save the uploaded file
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            st.success(f"File {uploaded_file.name} uploaded successfully!")

            # Trigger Seurat R script
            st.write("Running Seurat analysis...")
            try:
                os.makedirs("seurat_outputs", exist_ok=True)  # Ensure output directory exists

                # Define output RDS file path
                output_rds_file = f"seurat_outputs/{uploaded_file.name}_results.rds"

                # Run the R script with both arguments: input file and output file
                subprocess.run(["Rscript", "process_seurat.R", file_path, output_rds_file], check=True)

                st.success("Seurat analysis completed successfully!")

                # Display the UMAP plot
                st.image("seurat_outputs/umap_singleR_plot.png", caption="UMAP Plot", use_container_width=True)
            except Exception as e:
                st.error(f"Error running Seurat: {e}")

elif menu == "Results":
    st.header("Results")
    st.write("View results from RNAseq or scRNAseq analysis pipelines here.")
    # Dynamically check for results and display plots side by side
    plot_paths = [
        "seurat_outputs/umap_singleR_plot.png",
        "seurat_outputs/umap_plot.png",
        "seurat_outputs/dot_plot.png",
        "seurat_outputs/violin_plot_HBA1.png",
        "seurat_outputs/violin_plot_LYZ.png",
        "seurat_outputs/violin_plot_PPBP.png", 
        "seurat_outputs/violin_plot_S100A8.png", 
        "seurat_outputs/violin_plot_S100A9.png",
        "seurat_outputs/violin_plot_S100A12.png", 
        "seurat_outputs/feature_plot_HBA1.png",
        "seurat_outputs/feature_plot_LYZ.png",
        "seurat_outputs/feature_plot_PPBP.png", 
        "seurat_outputs/feature_plot_S100A8.png", 
        "seurat_outputs/feature_plot_S100A9.png",
        "seurat_outputs/feature_plot_S100A12.png",# Example for GeneA
        "seurat_outputs/cluster_heatmap.png",
        "seurat_outputs/pca_elbow_plot.png"  
    ]

    available_plots = [path for path in plot_paths if os.path.exists(path)]

    if available_plots:
        st.write("### Generated Plots:")

        # Create columns dynamically based on the number of plots
        num_cols = min(len(available_plots), 3)  # Show up to 3 plots in a row
        cols = st.columns(num_cols)

        for i, plot_path in enumerate(available_plots):
            col = cols[i % num_cols]  # Rotate through the columns
            with col:
                st.image(plot_path, caption=os.path.basename(plot_path), use_container_width=True)
    else:
        st.warning("No results available. Please upload and analyze data first.")


elif menu == "Analysis":
    st.header("Data Analysis")

    # Dynamically check for Seurat outputs in the 'seurat_outputs' folder
    results_folder = "seurat_outputs"

    # Check for pie chart or create it dynamically
    pie_chart_path = os.path.join(results_folder, "cell_type_statistics.png")
    cell_type_summary_file = os.path.join(results_folder, "cell_type_summary.csv")
    
    if os.path.exists(cell_type_summary_file):
        # Load cell type summary
        st.subheader("Cell-Type Statistics")
        cell_type_data = pd.read_csv(cell_type_summary_file)

        # Generate pie chart dynamically
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.pie(
            cell_type_data['Count'], 
            labels=cell_type_data['Cell Type'], 
            autopct='%1.1f%%', 
            startangle=140
        )
        ax.set_title("Cell-Type Distribution")
        st.pyplot(fig)

        # Display summary table
        st.write("### Cell-Type Summary")
        st.dataframe(cell_type_data)
    else:
        st.warning("No cell-type summary available. Please ensure analysis is complete.")

    # Check for differential genes data
    differential_genes_file = os.path.join(results_folder, "differential_genes.csv")
    if os.path.exists(differential_genes_file):
        st.subheader("Top Differential Genes for Each Cluster")
        diff_genes_data = pd.read_csv(differential_genes_file)

        # Display table
        st.dataframe(diff_genes_data)

        # Add search/filter functionality
        search_term = st.text_input("Search Genes or Clusters")
        if search_term:
            filtered_data = diff_genes_data[
                diff_genes_data.apply(
                    lambda row: search_term.lower() in row.to_string().lower(), axis=1
                )
            ]
            st.write("Filtered Results:")
            st.dataframe(filtered_data)
    else:
        st.warning("No differential genes data available.")

    # Insights Section
    st.write("### Insights")
    st.write("""
    - **Cell-Type Pie Chart**: Provides a visual representation of the proportion of cell types identified in the dataset.
    - **Top Differential Genes**: Highlights the most differentially expressed genes for each cluster.
    - **UMAP Plot**: Allows visualization of cell clusters and their annotations.
    """)

    # Display additional plots from Seurat outputs (e.g., UMAP, heatmaps, feature plots)
    available_plots = [f for f in os.listdir(results_folder) if f.endswith(".png") and "umap" in f]
    if available_plots:
        st.write("### Additional Plots")
        for plot_file in available_plots:
            st.image(os.path.join(results_folder, plot_file), caption=plot_file, use_container_width=True)
    else:
        st.warning("No additional plots available.")
    
elif menu == "Neo Dashboard":
    st.header("Neocellomics Dashboard")
    st.write("""
        Welcome to the Neocellomics Dashboard! This interactive dashboard provides advanced visualizations for patient profiles, single-cell results, and transcriptomics analysis.
    """)

    # Embed the Dash app
    DASH_URL = "http://127.0.0.1:8050"  # Replace with your hosted Dash app URL if deployed remotely
    st.components.v1.iframe(src=DASH_URL, height=800, scrolling=True)


elif menu == "Reports":
    st.header("Generate Reports")
    results_folder = "seurat_outputs"

    # Gather all plots and required data
    available_plots = [f for f in os.listdir(results_folder) if f.endswith(".png")]
    cell_type_summary_file = os.path.join(results_folder, "cell_type_summary.csv")
    differential_genes_file = os.path.join(results_folder, "differential_genes.csv")

    # Display plots for inclusion in the report
    if available_plots:
        st.write("### Plots to Include in Report:")
        for plot_file in available_plots:
            st.image(os.path.join(results_folder, plot_file), caption=plot_file, use_container_width=True)

    # Load cell type summary and differential genes data
    cell_type_data = None
    if os.path.exists(cell_type_summary_file):
        cell_type_data = pd.read_csv(cell_type_summary_file)
        st.write("### Cell-Type Summary Table:")
        st.dataframe(cell_type_data)

    diff_genes_data = None
    if os.path.exists(differential_genes_file):
        diff_genes_data = pd.read_csv(differential_genes_file)
        st.write("### Top Differential Genes Table:")
        st.dataframe(diff_genes_data)

    # Generate PDF Report
    if st.button("Generate PDF Report"):
        from fpdf import FPDF
        import matplotlib.pyplot as plt

        # Initialize PDF
        pdf = FPDF()
        pdf.set_auto_page_break(auto=True, margin=15)

        # Add Cover Page
        pdf.add_page()
        pdf.set_font("Arial", size=16)
        pdf.cell(200, 10, txt="CellGuide Transcriptomic Report", ln=True, align="C")
        pdf.ln(20)

        # Add Cell-Type Summary Pie Chart
        if cell_type_data is not None:
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.pie(
                cell_type_data['Count'],
                labels=cell_type_data['Cell Type'],
                autopct='%1.1f%%',
                startangle=140
            )
            ax.set_title("Cell-Type Distribution")
            pie_chart_path = os.path.join(results_folder, "cell_type_pie_chart_temp.png")
            fig.savefig(pie_chart_path)
            plt.close(fig)

            pdf.add_page()
            pdf.set_font("Arial", size=14)
            pdf.cell(200, 10, txt="Cell-Type Distribution", ln=True, align="L")
            pdf.image(pie_chart_path, x=10, y=40, w=180)

        # Add Cell-Type Summary Table
        if cell_type_data is not None:
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            pdf.cell(200, 10, txt="Cell-Type Summary Table", ln=True, align="L")
            pdf.ln(10)
            for _, row in cell_type_data.iterrows():
                pdf.cell(200, 10, txt=f"{row['Cell Type']}: {row['Count']} cells", ln=True, align="L")

        # Add Differential Genes Table
        if diff_genes_data is not None:
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            pdf.cell(200, 10, txt="Top Differential Genes Table", ln=True, align="L")
            pdf.ln(10)
            for i, row in diff_genes_data.iterrows():
                if i >= 20:  # Limit rows for brevity
                    pdf.cell(200, 10, txt="... (truncated)", ln=True, align="L")
                    break
                pdf.cell(200, 10, txt=f"Cluster {row['cluster']}, Gene: {row['gene']}, Log2FC: {row['avg_log2FC']}", ln=True, align="L")

        # Add Plots
        for plot_file in available_plots:
            plot_path = os.path.join(results_folder, plot_file)
            pdf.add_page()
            pdf.set_font("Arial", size=12)
            pdf.cell(200, 10, txt=f"Plot: {plot_file}", ln=True, align="L")
            pdf.image(plot_path, 
            x=10, y=40, w=180)

        # Save the PDF in memory
        pdf_output_path = "CellGuide_Transcriptomic_Report.pdf"
        pdf.output(pdf_output_path)

        # Read the PDF file as binary
        with open(pdf_output_path, "rb") as pdf_file:
            pdf_data = pdf_file.read()

        # Provide download link for the PDF
        st.success("PDF Report generated successfully!")
        st.download_button(
            label="Download Report",
            data=pdf_data,
            file_name="CellGuide_Transcriptomic_Report.pdf",
            mime="application/pdf"
        )
    else:
        st.warning("No data available to generate a report.")


