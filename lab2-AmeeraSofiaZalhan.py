import pandas as pd
import networkx as nx
import streamlit as st
import matplotlib.pyplot as plt
import requests

# Function to retrieve PPI data from BioGRID
def retrieve_ppi_biogrid(target_protein):
    api_key = '571ad712ab813d1eb1f93d9ade0bdda6'
    url = f"https://webservice.thebiogrid.org/interactions/?searchNames=true&geneList={target_protein}&format=json&accesskey={api_key}"
    response = requests.get(url)
    data = response.json()
    if data:
        df = pd.DataFrame(data).T  # Transpose if data is nested by interaction ID
        if 'OFFICIAL_SYMBOL_A' in df.columns and 'OFFICIAL_SYMBOL_B' in df.columns:
            df = df.rename(columns={'OFFICIAL_SYMBOL_A': 'protein1', 'OFFICIAL_SYMBOL_B': 'protein2'})
        return df[['protein1', 'protein2']]
    else:
        return pd.DataFrame(columns=['protein1', 'protein2'])

# Function to retrieve PPI data from STRING
def retrieve_ppi_string(target_protein):
    url = f"https://string-db.org/api/json/network?identifiers={target_protein}&species=9606"
    response = requests.get(url)
    data = response.json()
    if data:
        df = pd.DataFrame(data)
        if 'preferredName_A' in df.columns and 'preferredName_B' in df.columns:
            df = df.rename(columns={'preferredName_A': 'protein1', 'preferredName_B': 'protein2'})
        return df[['protein1', 'protein2']]
    else:
        return pd.DataFrame(columns=['protein1', 'protein2'])

# Function to generate a network from PPI DataFrame
def generate_network(df):
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['protein1'], row['protein2'])
    return G

# Function to calculate centrality measures
def get_centralities(G):
    centralities = {
        "Degree": nx.degree_centrality(G),
        "Betweenness": nx.betweenness_centrality(G),
        "Closeness": nx.closeness_centrality(G),
        "PageRank": nx.pagerank(G)
    }
    
    # handling convergence issues
    try:
        centralities["Eigenvector"] = nx.eigenvector_centrality(G, max_iter=500, tol=1e-06)
    except nx.PowerIterationFailedConvergence:
        centralities["Eigenvector"] = "Convergence failed - Unable to compute"

    return centralities


# Streamlit app
def main():
    st.set_page_config(page_title="Lab 2 - PPI Network Analysis", layout="wide")
    st.title("Lab 2 - Protein-Protein Interaction (PPI)")
    
    # Protein input and database selection on main page
    st.subheader("🔎 Input Parameters")
    target_protein = st.text_input("Enter Protein ID:")
    database_choice = st.selectbox("Choose Database:", ["BioGRID", "STRING"])

    if st.button("Retrieve PPI Data"):
        if database_choice == "BioGRID":
            df = retrieve_ppi_biogrid(target_protein)
        else:
            df = retrieve_ppi_string(target_protein)
        
        if not df.empty:
            # Limiting to the first 10 rows for display
            df_display = df.head(10)
            col1, col2 = st.columns(2, gap="large")

            # Column 1: PPI data and network visualization
            with col1:
                st.subheader("🧬 PPI Data Information")
                with st.expander("View PPI Data (First 10 rows)", expanded=True):
                    st.dataframe(df_display)

                st.write("**Network Statistics**")
                st.metric(label="Edges in Displayed Data", value=len(df_display))
                st.metric(label="Nodes in Displayed Data", value=len(set(df_display['protein1']).union(df_display['protein2'])))
                
                # Generate and display the network graph
                G = generate_network(df_display)
                fig, ax = plt.subplots()
                nx.draw(G, with_labels=True, node_color="skyblue", edge_color="gray", node_size=500, font_size=10, ax=ax)
                st.pyplot(fig)

            # Column 2: Centrality Measures
            with col2:
                st.subheader("📊 Centrality Measures")
                centralities = get_centralities(G)
                
                for measure, values in centralities.items():
                    with st.expander(f"{measure.capitalize()} Centrality", expanded=True):
                        if isinstance(values, dict):
                            st.json(values)
                        else:
                            st.write(values)

# Simple footer using st.markdown
    st.markdown("---")
    st.markdown("Lab 2 - Ameera Sofia Zalhan")

if __name__ == "__main__":
    main()
