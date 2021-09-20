def add_find_candidates_arguments(parser):
    """
    Add arguments to a parser for sub-command "stitch"
    :param parser: argeparse object
    :return:
    """
    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        required=True,
        help="Path to directory containing HDF files."
    )
    parser.add_argument(
        "-b",
        "--bam",
        type=str,
        required=True,
        help="BAM file containing mapping between reads and the draft assembly."
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="Input reference/assembly file."
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        type=str,
        required=True,
        help="Name of the sample."
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=str,
        required=True,
        help="Path to output directory."
    )
    parser.add_argument(
        "-t",
        "--threads",
        required=True,
        type=int,
        help="Number of threads."
    )
    parser.add_argument(
        "-hp",
        "--use_hp_info",
        default=False,
        action='store_true',
        help="If set then haplotype-aware mode will be enabled."
    )
    parser.add_argument(
        "--allowed_multiallelics",
        type=int,
        required=False,
        default=4,
        help="Number of maximum multialleleic variants allowed per site. Default is 4"
    )
    parser.add_argument(
        "--snp_p_value",
        required=False,
        type=float,
        default=0.4,
        help="Predicted value used for a SNP to be considered a candidate. Default is 0.4"
    )
    parser.add_argument(
        "--insert_p_value",
        required=False,
        type=float,
        default=0.1,
        help="Predicted value used for a insert to be considered a candidate. Default is 0.1"
    )
    parser.add_argument(
        "--delete_p_value",
        required=False,
        type=float,
        default=0.2,
        help="Predicted value used for a delete to be considered a candidate. Default is 0.2"
    )
    parser.add_argument(
        "--freq_based",
        default=False,
        action='store_true',
        help="If set then frequency based variants will be reported."
    )
    parser.add_argument(
        "--freq",
        required=False,
        type=float,
        default=0.10,
        help="If frequency based variant finding in enabled then this frequency will be the threshold."
    )
    return parser