import streamlit as st
import bio_algorithms as bio

# ================== PAGE CONFIG ==================
st.set_page_config(
    page_title="Bioinformatics Toolkit",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ================== CUSTOM STYLING ==================
st.markdown("""
<style>
    /* Main container */
    .main .block-container {
        padding-top: 2rem;
        padding-bottom: 2rem;
        max-width: 1200px;
    }
    
    /* Header styling */
    .stApp h1 {
        background: linear-gradient(120deg, #1e88e5 0%, #7c4dff 50%, #00bfa5 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        font-size: 2.5rem !important;
        font-weight: 700 !important;
        margin-bottom: 1rem;
    }
    
    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #0d1117 0%, #161b22 100%);
    }
    
    [data-testid="stSidebar"] .stRadio label {
        color: #c9d1d9 !important;
        padding: 0.5rem 1rem;
        border-radius: 0.5rem;
        transition: all 0.3s ease;
    }
    
    [data-testid="stSidebar"] .stRadio label:hover {
        background: rgba(99, 102, 241, 0.1);
    }
    
    /* Result cards */
    .result-card {
        background: linear-gradient(135deg, #1a1f35 0%, #0d1117 100%);
        border: 1px solid #30363d;
        border-radius: 1rem;
        padding: 1.5rem;
        margin: 1rem 0;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
    }
    
    .result-value {
        font-size: 2rem;
        font-weight: 700;
        background: linear-gradient(90deg, #00d4ff, #7c4dff);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    
    /* Sequence display */
    .sequence-box {
        background: #0d1117;
        border: 1px solid #30363d;
        border-radius: 0.75rem;
        padding: 1rem;
        font-family: 'Monaco', 'Menlo', monospace;
        font-size: 0.9rem;
        color: #58a6ff;
        word-wrap: break-word;
        overflow-x: auto;
        margin: 0.5rem 0;
    }
    
    /* Buttons */
    .stButton > button {
        background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 100%);
        color: white;
        border: none;
        border-radius: 0.75rem;
        padding: 0.75rem 2rem;
        font-weight: 600;
        font-size: 1rem;
        transition: all 0.3s ease;
        box-shadow: 0 4px 15px rgba(99, 102, 241, 0.4);
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(99, 102, 241, 0.6);
    }
    
    /* Text areas */
    .stTextArea textarea {
        background: #0d1117;
        border: 2px solid #30363d;
        border-radius: 0.75rem;
        color: #c9d1d9;
        font-family: 'Monaco', 'Menlo', monospace;
    }
    
    .stTextArea textarea:focus {
        border-color: #6366f1;
        box-shadow: 0 0 0 3px rgba(99, 102, 241, 0.2);
    }
    
    /* File uploader */
    [data-testid="stFileUploader"] {
        background: #0d1117;
        border: 2px dashed #30363d;
        border-radius: 0.75rem;
        padding: 1rem;
    }
    
    /* Info/Success/Error boxes */
    .stAlert {
        border-radius: 0.75rem;
    }
    
    /* Tabs styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 0.5rem;
    }
    
    .stTabs [data-baseweb="tab"] {
        background: transparent;
        border-radius: 0.5rem;
        padding: 0.5rem 1rem;
        transition: all 0.3s ease;
    }
    
    /* Divider */
    hr {
        border: none;
        height: 1px;
        background: linear-gradient(90deg, transparent, #30363d, transparent);
        margin: 2rem 0;
    }
    
    /* Hero section */
    .hero-badge {
        display: inline-block;
        background: linear-gradient(135deg, rgba(99, 102, 241, 0.2), rgba(139, 92, 246, 0.2));
        border: 1px solid rgba(99, 102, 241, 0.3);
        border-radius: 2rem;
        padding: 0.5rem 1rem;
        font-size: 0.85rem;
        color: #a5b4fc;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)


# ================== HELPER FUNCTIONS ==================
def get_sequence_from_input(text_input, uploaded_file, input_type):
    """Return sequence from text input or FASTA file"""
    if input_type == "Type Sequence":
        seq = text_input.strip().upper()
        if not seq:
            raise ValueError("Please enter a sequence")
        return seq
    else:
        if uploaded_file is None:
            raise ValueError("Please upload a FASTA file")
        content = uploaded_file.read().decode("utf-8")
        lines = content.strip().split("\n")
        seq_lines = [line.strip() for line in lines if not line.startswith(">")]
        return "".join(seq_lines).upper()


def display_result(label, value, icon=""):
    """Display result in a styled card"""
    st.markdown(f"""
    <div class="result-card">
        <p style="color: #8b949e; margin-bottom: 0.5rem;">{label}</p>
        <p class="result-value">{value}</p>
    </div>
    """, unsafe_allow_html=True)


def display_sequence(label, sequence, icon="🧬"):
    """Display sequence in a styled box"""
    st.markdown(f"**{icon} {label}:**")
    st.markdown(f'<div class="sequence-box">{sequence}</div>', unsafe_allow_html=True)


def create_input_block(key_prefix, label="Sequence", show_file_upload=True):
    """Create standardized input block with text area and optional file upload"""
    input_type = st.radio(
        f"Input method for {label}:",
        ["Type Sequence", "📁 Upload FASTA"],
        key=f"{key_prefix}_input_type",
        horizontal=True
    )
    
    text_input = ""
    uploaded_file = None
    
    if input_type == "Type Sequence":
        text_input = st.text_area(
            f"Enter {label}:",
            height=120,
            key=f"{key_prefix}_text",
            placeholder="Enter DNA/RNA sequence (e.g., ATGCATGCATGC)"
        )
    else:
        uploaded_file = st.file_uploader(
            f"Upload FASTA file for {label}:",
            type=["fa", "fasta"],
            key=f"{key_prefix}_file"
        )
    
    return text_input, uploaded_file, input_type


# ================== SIDEBAR NAVIGATION ==================
with st.sidebar:
    st.markdown("## Navigation")
    st.markdown("---")
    
    tool = st.radio(
        "Select Tool:",
        [
            "Home",
            "GC Content",
            "Reverse / Complement",
            "DNA Translation",
            "Hamming Distance",
            "Edit Distance",
            "Suffix Array",
            "Overlap Graph",
            "Boyer-Moore Search"
        ],
        label_visibility="collapsed"
    )
    
    # st.markdown("---")
    # st.markdown("""
    # <div style="text-align: center; color: #8b949e; font-size: 0.8rem;">
    #     <p>Bioinformatics Toolkit</p>
    # </div>
    # """, unsafe_allow_html=True)


# ================== HOME PAGE ==================
if tool == "Home":
    # st.markdown('<span class="hero-badge">Bioinformatics Analysis Suite</span>', unsafe_allow_html=True)
    st.title("Welcome to the Bioinformatics Toolkit")
    
    st.markdown("""
    A comprehensive toolkit for DNA/RNA sequence analysis featuring powerful algorithms 
    for molecular biology research and education.
    """)
    
    st.markdown("---")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.markdown("""
        #### Sequence Analysis
        - GC Content
        - Reverse/Complement
        - DNA Translation
        """)
    
    with col2:
        st.markdown("""
        #### Distance Metrics
        - Hamming Distance
        - Edit Distance
        """)
    
    with col3:
        st.markdown("""
        #### Data Structures
        - Suffix Array
        - Overlap Graph
        """)
    
    with col4:
        st.markdown("""
        #### Pattern Search
        - Boyer-Moore Algorithm
        """)
# ================== GC CONTENT ==================
elif tool == "GC Content":
    st.title("GC Content Calculator")
    st.markdown("Calculate the proportion of guanine (G) and cytosine (C) bases in your DNA sequence.")
    
    st.markdown("---")
    
    text_input, uploaded_file, input_type = create_input_block("gc")
    
    if st.button("Calculate GC Content", type="primary", use_container_width=True):
        try:
            seq = get_sequence_from_input(text_input, uploaded_file, input_type)
            result = bio.GC_Content(seq)
            display_result("GC Content", f"{result:.4f} ({result*100:.2f}%)")
            st.success(f"Analyzed sequence of length {len(seq)} bp")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== REVERSE / COMPLEMENT ==================
elif tool == "Reverse / Complement":
    st.title("Reverse / Complement")
    st.markdown("Generate reverse, complement, and reverse complement of your DNA sequence.")
    
    st.markdown("---")
    
    text_input, uploaded_file, input_type = create_input_block("rc")
    
    if st.button("Generate All", type="primary", use_container_width=True):
        try:
            seq = get_sequence_from_input(text_input, uploaded_file, input_type)
            
            col1, col2 = st.columns(2)
            
            with col1:
                display_sequence("Original", seq)
                display_sequence("Reverse", bio.Reverse(seq))
            
            with col2:
                display_sequence("Complement", bio.Complement(seq))
                display_sequence("Reverse Complement", bio.Reverse_Complement(seq))
            
            st.success(f"Processed sequence of length {len(seq)} bp")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== DNA TRANSLATION ==================
elif tool == "DNA Translation":
    st.title("DNA to Protein Translation")
    st.markdown("Translate your DNA sequence into amino acids using the standard genetic code.")
    
    st.markdown("---")
    
    text_input, uploaded_file, input_type = create_input_block("tr")
    
    if st.button("Translate", type="primary", use_container_width=True):
        try:
            seq = get_sequence_from_input(text_input, uploaded_file, input_type)
            protein = bio.translate_dna(seq)
            
            display_sequence("DNA Sequence", seq)
            display_sequence("Protein Sequence", protein)
            
            st.info(f"Translated {len(seq)} nucleotides → {len(protein)} amino acids")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== HAMMING DISTANCE ==================
elif tool == "Hamming Distance":
    st.title("Hamming Distance")
    st.markdown("Calculate the number of positions where corresponding symbols differ between two equal-length sequences.")
    
    st.markdown("---")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Sequence 1")
        text1, file1, type1 = create_input_block("ham1", "Sequence 1")
    
    with col2:
        st.markdown("### Sequence 2")
        text2, file2, type2 = create_input_block("ham2", "Sequence 2")
    
    if st.button("Calculate Distance", type="primary", use_container_width=True):
        try:
            seq1 = get_sequence_from_input(text1, file1, type1)
            seq2 = get_sequence_from_input(text2, file2, type2)
            distance = bio.hamming_distance(seq1, seq2)
            
            display_result("Hamming Distance", distance)
            st.success(f"Compared sequences of length {len(seq1)} bp")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== EDIT DISTANCE ==================
elif tool == "Edit Distance":
    st.title("Edit Distance (Levenshtein)")
    st.markdown("Calculate the minimum number of single-character edits required to transform one sequence into another.")
    
    st.markdown("---")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### Sequence 1")
        text1, file1, type1 = create_input_block("ed1", "Sequence 1")
    
    with col2:
        st.markdown("### Sequence 2")
        text2, file2, type2 = create_input_block("ed2", "Sequence 2")
    
    if st.button("Calculate Edit Distance", type="primary", use_container_width=True):
        try:
            seq1 = get_sequence_from_input(text1, file1, type1)
            seq2 = get_sequence_from_input(text2, file2, type2)
            distance = bio.edit_distance(seq1, seq2)
            
            display_result("Edit Distance", distance)
            st.info(f"Comparing sequences: {len(seq1)} bp vs {len(seq2)} bp")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== SUFFIX ARRAY ==================
elif tool == "Suffix Array":
    st.title("Suffix Array")
    st.markdown("Construct a suffix array and its inverse for efficient string operations.")
    
    st.markdown("---")
    
    text_input, uploaded_file, input_type = create_input_block("sa")
    
    if st.button("Build Suffix Array", type="primary", use_container_width=True):
        try:
            seq = get_sequence_from_input(text_input, uploaded_file, input_type)
            sa = bio.suffix_array(seq)
            inv = bio.inverse_suffix_array(sa)
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("### Suffix Array")
                st.code(str(sa), language=None)
            
            with col2:
                st.markdown("### Inverse Suffix Array")
                st.code(str(inv), language=None)
            
            st.success(f"Built suffix array for sequence of length {len(seq)}")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== OVERLAP GRAPH ==================
elif tool == "Overlap Graph":
    st.title("Overlap Graph")
    st.markdown("Build an overlap graph from a set of reads for sequence assembly.")
    
    st.markdown("---")
    
    input_type = st.radio(
        "Input method:",
        ["Type Reads", "📁 Upload FASTA"],
        key="og_input_type",
        horizontal=True
    )
    
    reads_input = []
    
    if input_type == "Type Reads":
        reads_text = st.text_area(
            "Enter comma-separated reads:",
            height=120,
            key="og_text",
            placeholder="ATGCG, GCGTC, GTCAA, CAAATG"
        )
    else:
        uploaded_file = st.file_uploader(
            "Upload FASTA file with reads:",
            type=["fa", "fasta"],
            key="og_file"
        )
    
    k = st.number_input("Minimum overlap (k):", min_value=1, value=3, step=1)
    
    if st.button("Build Graph", type="primary", use_container_width=True):
        try:
            if input_type == "Type Reads":
                reads = [r.strip().upper() for r in reads_text.split(',') if r.strip()]
            else:
                if uploaded_file is None:
                    raise ValueError("Please upload a FASTA file")
                content = uploaded_file.read().decode("utf-8")
                records = []
                current_seq = ""
                for line in content.strip().split("\n"):
                    if line.startswith(">"):
                        if current_seq:
                            records.append(current_seq)
                        current_seq = ""
                    else:
                        current_seq += line.strip()
                if current_seq:
                    records.append(current_seq)
                reads = [seq.upper() for seq in records]
            
            if len(reads) < 2:
                raise ValueError("Please provide at least 2 reads")
            
            graph = bio.overlap_graph(reads, k)
            
            st.markdown("### Overlap Graph Result")
            
            if graph:
                for source, targets in graph.items():
                    with st.expander(f"{source[:30]}..." if len(source) > 30 else f"{source}"):
                        for target, overlap_len in targets:
                            st.markdown(f"→ **{target}** (overlap: {overlap_len})")
            else:
                st.warning("No overlaps found with the given minimum overlap length.")
            
            st.success(f"Analyzed {len(reads)} reads with minimum overlap k={k}")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== BOYER-MOORE SEARCH ==================
elif tool == "Boyer-Moore Search":
    st.title("Boyer-Moore Pattern Search")
    st.markdown("Find all occurrences of a pattern in a text using the efficient Boyer-Moore algorithm.")
    
    st.markdown("---")
    
    st.markdown("### Text (Haystack)")
    text_input, text_file, text_type = create_input_block("bm_text", "Text")
    
    st.markdown("### Pattern (Needle)")
    pattern_input, pattern_file, pattern_type = create_input_block("bm_pattern", "Pattern")
    
    if st.button("Search Pattern", type="primary", use_container_width=True):
        try:
            text = get_sequence_from_input(text_input, text_file, text_type)
            pattern = get_sequence_from_input(pattern_input, pattern_file, pattern_type)
            
            matches = bio.boyer_moore_search(text, pattern)
            
            if matches:
                display_result("Matches Found", len(matches))
                
                st.markdown("### Match Positions")
                st.markdown(f'<div class="sequence-box">{", ".join(map(str, matches))}</div>', 
                           unsafe_allow_html=True)
                
                # Show context around first few matches
                if len(matches) <= 5:
                    st.markdown("### Match Context")
                    for pos in matches:
                        start = max(0, pos - 10)
                        end = min(len(text), pos + len(pattern) + 10)
                        context = text[start:pos] + f"**{text[pos:pos+len(pattern)]}**" + text[pos+len(pattern):end]
                        st.markdown(f"Position {pos}: ...{context}...")
            else:
                st.warning("No matches found for the given pattern.")
            
            st.success(f"Searched pattern of length {len(pattern)} in text of length {len(text)}")
        except Exception as e:
            st.error(f"Error: {str(e)}")


# ================== FOOTER ==================
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #8b949e; padding: 1rem;">
    <p>🧬 Bioinformatics Toolkit</p>
</div>
""", unsafe_allow_html=True)
