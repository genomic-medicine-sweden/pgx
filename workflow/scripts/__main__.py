import click

@click.group(help="CLI tool for pgx")
def cli():
    print("HELP text")

@cli.command(short_help="get target variants")
def get_target_variants():
    print("get_target_variants.py")

@cli.command(short_help="append_rsid to gdf")
def get_append_rsid_to_gdf():
    print("append_rsid_to_gdf")

@click.option("-d", "--description", prompt=True, required=True, type=str, help="A short description of your pipeline")   
def get_append_rsid_to_gdf(description):
    input_gdf = snakemake.input["gdf"]
    target_bed = snakemake.params["target_bed"]
    output_file = snakemake.output["out"]   
    c_gdf = GDF(input_gdf)
    c_gdf.rsid_per_position(target_bed)
    c_gdf.write_proccessed_gdf(output_file)

if __name__ == '__main__':
    cli()

