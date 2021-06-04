# 导入 -> 系统 package
use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Statistics::R;
use Statistics::Descriptive;
use Excel::Writer::XLSX;
$|=1;
 
# 定义 -> 常量
use constant SCRIPTDIR => (File::Spec->splitpath(File::Spec->rel2abs($0)))[1];
use constant PWD => $ENV{"PWD"};

# 定义 -> 核心变量
 


# 检测 -> 脚本输入
my ($input, $motif_file, $format, $output_dir, $min_reads, $min_rc, $output_model, $if_help);
GetOptions(
    "input|i=s"         => \$input,
    "motif|m=s"         => \$motif_file,
    "format|f=s"        => \$format,
    "output_dir|o=s"    => \$output_dir,
    "min_reads=s"       => \$min_reads,
    "min_rc=s"          => \$min_rc,
    "output_model=s"    => \$output_model,
    "help|h"            => \$if_help,
);
die help() if(defined $if_help or (not defined $input or not defined $motif_file or not defined $format or not defined $output_dir));
die "[Error] format must be 'row' or 'col'" if($format ne 'row' and $format ne 'col');
###################################################################### 主程序
# 分型系统条件定义
$min_reads = 30 if(not defined $min_reads or $min_reads !~ /^\d+$/);
$min_rc    = 0 if(not defined $min_rc     or $min_rc !~ /^\d+\.\d+$/);
$output_model = 'normal' if(not defined $output_model);
die "[Error] output_model can only be normal/complete/all\n" if($output_model ne 'normal' and $output_model ne 'complete' and $output_model ne 'all'); 

# 输出目录创建
mkdir $output_dir if(not -e $output_dir);

########
# （1）读取reads数量与频率
########
my @samples;  # 样本列表
my %hashSTR;  # reads数、频率
my %hashMotif;  # 参考片段的motif信息
read_motif($motif_file, \%hashMotif);
read_str_row($input, \%hashSTR, \%hashMotif, \@samples, $min_reads) if($format  eq 'row'); 
read_str_col($input, \%hashSTR, \%hashMotif, \@samples, $min_reads) if($format  eq 'col'); 

########
# （2）前滑移校正
########
my %hashSlipRatio = calculate_slip_ratio(\%hashSTR, \%hashMotif, \@samples); # 计算滑移比例
my %hashSTRClean  = slip_correction(\%hashSTR, \%hashMotif, \%hashSlipRatio, \@samples); # 校正
 
########
# （3）扩增校正比例
######## 
my %hashAMPR = calculate_amplification_ratio(\%hashSTRClean, \%hashMotif, \@samples, $min_rc); # 计算扩增比例 

########
# （4）分型
########   
my %hashTyping = typing(\%hashSTRClean, \%hashMotif, \%hashAMPR, \@samples, $min_rc); # 分型
    
########
# （5）结果输出
########  
print "output result ... ";
# 5.1 TXT 文件生成
my $str_reads_raw      = "$output_dir/1.str.reads.raw.txt";           # 1. 原始reads数
my $str_freq_raw       = "$output_dir/2.str.freq.raw.txt";            # 2. 原始freq
my $str_slip_ratio     = "$output_dir/3.str.slip.ratio.txt";          # 3. 滑移比例
my $str_reads_correct  = "$output_dir/4.str.reads.correct.txt";       # 4. 滑移校正后的reads数
my $str_freq_correct   = "$output_dir/5.str.freq.correct.txt";        # 5. 滑移校正后的频率
my $str_ampr           = "$output_dir/6.str.amplification.ratio.txt"; # 6. 扩增校正比例
my $str_type_relative  = "$output_dir/7.str.typing.relative.txt";     # 7. 分型相对拷贝数
my $str_type_absolute  = "$output_dir/8.str.typing.absolute.txt";     # 8. 分型绝对拷贝数
my $str_genotype       = "$output_dir/9.str.genotype.txt";            # 9. 分型结果
my $str_allele         = "$output_dir/10.str.allele.txt";             # 10. 等位基因
my $str_allele_r       = "$output_dir/11.str.allele.relative.txt";    # 11. 等位基因相对拷贝数
my $log_file           = "$output_dir/log.txt";                       # 输出信息的log记录
 

open LOG, ">>$log_file";
output_reads_freq(\%hashSTR, \%hashMotif, 'Reads', \@samples, $str_reads_raw);          # 1
output_reads_freq(\%hashSTR, \%hashMotif, 'Freq', \@samples, $str_freq_raw);            # 2
output_slip_ratio(\%hashSlipRatio, \%hashMotif, $str_slip_ratio);                                    # 3
output_reads_freq(\%hashSTRClean, \%hashMotif, 'Reads', \@samples, $str_reads_correct); # 4
output_reads_freq(\%hashSTRClean, \%hashMotif, 'Freq', \@samples, $str_freq_correct);   # 5
output_ampr(\%hashAMPR, \%hashMotif, $str_ampr);                                                     # 6
output_type_relative_absolute(\%hashTyping, \%hashMotif, 'Relative', \@samples, $str_type_relative); # 7
output_type_relative_absolute(\%hashTyping, \%hashMotif, 'Absolute', \@samples, $str_type_absolute); # 8
output_genotype(\%hashTyping, \%hashMotif, 'genotype', \@samples, $str_genotype); # 9
output_genotype(\%hashTyping, \%hashMotif, 'allele', \@samples, $str_allele);     # 10 
output_genotype(\%hashTyping, \%hashMotif, 'allele_r', \@samples, $str_allele_r);     # 10 
output_excel();  # 5.2 TXT 转 excel(注意：由于涉及到的变量过多，传递参数较为麻烦，故使用全局变量的形式执行)
close LOG;

# 清除中间结果
print "OK\n\n";
print "result file: $output_dir/SSRSeq_type.xlsx\n";

###################################################################### 子函数 ###############################################################

sub output_excel{
    my $excel          = "$output_dir/SSRSeq_type.xlsx";          # 提供给客户
    my $excel_complete = "$output_dir/SSRSeq_type_complete.xlsx"; # 生信内部阅读,包含Ycontrol样本

    print LOG "output $excel\n" if($output_model eq 'all' or $output_model eq 'normal');
    print LOG "output $excel_complete\n" if($output_model eq 'all' or $output_model eq 'complete');

    print LOG "output process : \n";
    my $workbook          = Excel::Writer::XLSX->new($excel) if($output_model eq 'all' or $output_model eq 'normal');# 结果文件
    my $workbook_complete = Excel::Writer::XLSX->new($excel_complete) if($output_model eq 'all' or $output_model eq 'complete');# 结果文件

    my %format            = sheet_format($workbook) if($output_model eq 'all' or $output_model eq 'normal'); 
    my %format_complete   = sheet_format($workbook_complete) if($output_model eq 'all' or $output_model eq 'complete'); 

    # (1) 原始reads数量输出
    print LOG "    1 - raw reads\n";
    reads_to_excel($workbook,          \%format,          'Reads', $str_reads_raw, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    reads_to_excel($workbook_complete, \%format_complete, 'Reads', $str_reads_raw, \@samples) if($output_model eq 'all' or $output_model eq 'complete');
    # (2) 原始频率输出
    print LOG "    2 - raw freq\n";
    freq_to_excel($workbook,          \%format,          'Freq', $str_freq_raw, \%hashMotif, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    freq_to_excel($workbook_complete, \%format_complete, 'Freq', $str_freq_raw, \%hashMotif, \@samples) if($output_model eq 'all' or $output_model eq 'complete');

    # (3) 滑移校正数据
    print LOG "    3 - slip ratio\n";
    txt_to_excel($workbook_complete, \%format_complete, 'SlipRatio', $str_slip_ratio) if($output_model eq 'all' or $output_model eq 'complete');

    # (4) 校正reads数量输出
    print LOG "    4 - correct reads\n";
    # reads_to_excel($workbook, \%format, 'Reads Correct', $str_reads_correct, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    reads_to_excel($workbook_complete, \%format_complete, 'Reads Correct', $str_reads_correct, \@samples) if($output_model eq 'all' or $output_model eq 'complete');

    # (5) 校正freq数量输出
    print LOG "    5 - correct freq\n";
    freq_to_excel($workbook, \%format, 'Freq Correct', $str_freq_correct, \%hashMotif, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    freq_to_excel($workbook_complete, \%format_complete, 'Freq Correct', $str_freq_correct, \%hashMotif, \@samples) if($output_model eq 'all' or $output_model eq 'complete');
    
    # (6) 滑移校正数据
    print LOG "    6 - amplification ratio\n";
    txt_to_excel($workbook_complete, \%format_complete, 'Amplification Ratio', $str_ampr) if($output_model eq 'all' or $output_model eq 'complete');

    # (7) 相对拷贝数
    print LOG "    7 - typing relative\n";
    # type_relative_absolute_to_excel($workbook,          \%format,          'Typing Relative', $str_type_relative, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    txt_to_excel($workbook_complete,                    \%format_complete, 'Typing Relative', $str_type_relative) if($output_model eq 'all' or $output_model eq 'complete');

    # (8) 绝对拷贝数
    print LOG "    8 - typing absolute\n";
    # type_relative_absolute_to_excel($workbook,          \%format,          'Typing Absolute', $str_type_absolute, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    txt_to_excel($workbook_complete,                    \%format_complete, 'Typing Absolute', $str_type_absolute) if($output_model eq 'all' or $output_model eq 'complete');

    # (9) 分型
    print LOG "    9 - genotype\n";
    genotype_to_excel($workbook,          \%format,          'Genotype', $str_genotype, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    txt_to_excel($workbook_complete,      \%format_complete, 'Genotype', $str_genotype) if($output_model eq 'all' or $output_model eq 'complete');

    # (10) 等位基因
    print LOG "    10 - allele\n";
    genotype_to_excel($workbook,          \%format,          'Allele', $str_allele, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    txt_to_excel($workbook_complete,      \%format_complete, 'Allele', $str_allele) if($output_model eq 'all' or $output_model eq 'complete');

    # (11) 等位基因 相对量
    print LOG "    11 - allele relative\n";
    # genotype_to_excel($workbook,          \%format,          'Allele Relative', $str_allele_r, \@samples) if($output_model eq 'all' or $output_model eq 'normal');
    txt_to_excel($workbook_complete,      \%format_complete, 'Allele Relative', $str_allele_r) if($output_model eq 'all' or $output_model eq 'complete');

    # (12) readme
    print LOG "    11 - readme\n";
    readme($workbook,          \%format,          SCRIPTDIR.'/readme.txt') if($output_model eq 'all' or $output_model eq 'normal');
    readme($workbook_complete, \%format_complete, SCRIPTDIR.'/readme_complete.txt') if($output_model eq 'all' or $output_model eq 'complete');

    print LOG "    *Finished* \n";
}


# 输出readme
sub readme{
    my $workbook    = shift @_;
    my $format      = shift @_;
    my $readme_file = shift @_;
 
    my $sheet = $workbook->add_worksheet('ReadMe');

    $sheet->set_row(0, 65);
    $sheet->set_column('A:A', 35);
    $sheet->set_column('B:B', 25);
    $sheet->set_column('C:C', 110);
    my $row = 0;
    open FILE,$readme_file;
    while(<FILE>){
        $_=~ s/\^/\n/g;
        my @datas = split /\t/, $_;
        my $col = 0;
        foreach my $data(@datas)
        {
            my $text = Encode::decode("gb2312",$data);
            $sheet->write($row, $col, $text, $format->{'readme1'})    if($row == 0 and $col == 0);
            $sheet->write($row, $col, $text, $format->{'readme2'})    if($row == 0 and $col == 1);
            $sheet->write($row, $col, $text, $format->{'readme2tmp'}) if($row == 0 and $col == 2);

            $sheet->write($row, $col, $text, $format->{'readme3'})    if($row >  0 and $col == 0);
            $sheet->write($row, $col, $text, $format->{'readme4'})    if($row >  0 and $col == 1);
            $sheet->write($row, $col, $text, $format->{'readme5'})    if($row >  0 and $col == 2);
            $col++;
        }
        $row++;
    }
    close FILE;   
}


# 样本的分型结果放到excel中
sub genotype_to_excel{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $samples    = shift @_;

    my %hashData = read_file($input_file);

    my @heads = ('target', @$samples);

    my $sheet = $workbook->add_worksheet($sheet_name);
    $sheet->set_row(0, 60);
    $sheet->set_column(0, 0, 20);
 

    # 表头
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {   
        my @datas;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            push @datas, $data;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write_string($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;
    }    
}

# 相对、绝对分型结果输出
sub type_relative_absolute_to_excel{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $samples    = shift @_;

    my %hashData = read_file($input_file);
    my @heads    = split /\t/, $hashData{'HEAD'};

    my $sheet = $workbook->add_worksheet($sheet_name);
    $sheet->set_row(0, 60);
    $sheet->set_column(0, 0, 20);
  

    # 表头
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++; 

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {   
        my @datas;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            push @datas, $data;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;
    }   
}

# txt文档直接输出到excel
sub txt_to_excel{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;

    my %hashData = read_file($input_file);
    my @heads    = split /\t/, $hashData{'HEAD'};
    
    my $sheet = $workbook->add_worksheet($sheet_name);
    $sheet->set_row(0, 60);
    $sheet->set_column(0, 0, 20);
  

    # 表头
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++; 

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {   
        my @datas;
        my @colors;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            my $color = "normal";
               $color = "red"    if($head eq 'SlipRatio' and exists $hashData{'Value'}{$row_count}{'Type'} and $hashData{'Value'}{$row_count}{'Type'} eq 'unitary_quadratic_equation' );
            push @datas, $data;
            push @colors, $color;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write_string($row, $col, $datas[$col], $format->{$colors[$col]});
        }
        $row++;
    }
}

# freq输出
sub freq_to_excel{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $hashMotif  = shift @_;
    my $samples    = shift @_;
    my $ycontrols  = "";
       $ycontrols  = shift @_ if(exists $_[0]);

    my %hashData = read_file($input_file);

    my @heads = ('target', 'str', @$samples);
    push @heads, @$ycontrols if($ycontrols ne '');

    my $sheet = $workbook->add_worksheet($sheet_name);
    $sheet->set_row(0, 60);
    $sheet->set_column(0, 0, 20);
 

    # 表头
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {   
        my $target       = $hashData{'Value'}{$row_count}{'target'};
        my $noise_cutoff = $hashMotif->{$target}{'noise_cutoff'}; # 设定的噪音阈值
        my $type_cutoff  = $hashMotif->{$target}{'type_cutoff'};  # 设定的分型阈值

        my @datas;
        my @colors;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            my $color = "normal";
               $color = "orange" if($data =~ /\d/ and $head ne 'target' and $head ne "str" and $data >= $noise_cutoff );
               $color = "red"    if($data =~ /\d/ and $head ne 'target' and $head ne "str" and $data >= $type_cutoff );
            push @datas, $data;
            push @colors, $color;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{$colors[$col]});
        }
        $row++;
    }
}

# reads数量输出
sub reads_to_excel{
    my $workbook   = shift @_;
    my $format     = shift @_;
    my $sheet_name = shift @_;
    my $input_file = shift @_;
    my $samples    = shift @_;
    my $ycontrols  = "";
       $ycontrols  = shift @_ if(exists $_[0]);

    my %hashData = read_file($input_file);

    my @heads = ('target', 'str', @$samples);
    push @heads, @$ycontrols if($ycontrols ne '');

    my $sheet = $workbook->add_worksheet($sheet_name);
    $sheet->set_row(0, 60);
    $sheet->set_column(0, 0, 20);
 

    # 表头
    my $row = 0;
    foreach my $col(0..$#heads)
    {
        $sheet->write_string($row, $col, $heads[$col], $format->{'title'});
    }
    $row++;

    # 内容
    foreach my $row_count( sort {$a<=>$b} keys %{$hashData{'Value'}})
    {
        my @datas;
        foreach my $head(@heads)
        {
            my $data = exists $hashData{'Value'}{$row_count}{$head} ? $hashData{'Value'}{$row_count}{$head} : "";
            push @datas, $data;
        }
        foreach my $col(0..$#datas)
        {
            $sheet->write($row, $col, $datas[$col], $format->{'normal'});
        }
        $row++;
    }
}

# 读取文件
sub read_file{
    my $file = shift @_;
    my %hashData;

    my $row = 0;
    open INPUT, $file;
    my $line0 = <INPUT>;
       $line0 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line0;

    $row++;
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        my @datas = split /\t/, $_;
        foreach my $col(0..$#heads)
        {
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $hashData{'Value'}{$row}{$heads[$col]} = $value;
        }
        $row++;
    }
    close INPUT;

    $hashData{'HEAD'} = join "\t", @heads;
    return %hashData;
}

# 表格格式模版
sub sheet_format{
    my ($workbook)=@_;
    my %format=();
    $format{'title'} = $workbook->add_format();
    $format{'title'} ->set_align('center');
    $format{'title'} ->set_align('vcenter');
    $format{'title'} ->set_size(12);
    $format{'title'} ->set_font("Times New Roman");
    $format{'title'} ->set_border();
    $format{'title'} ->set_bg_color("yellow");
    $format{'title'} ->set_color("black");

    $format{'red'} = $workbook->add_format();
    $format{'red'} ->set_align('center');
    $format{'red'} ->set_align('vcenter');
    $format{'red'} ->set_size(12);
    $format{'red'} ->set_font("Times New Roman");
    $format{'red'} ->set_bg_color("#ff0000");
    $format{'red'} ->set_border();
    
    $format{'normal'} = $workbook->add_format();
    $format{'normal'} ->set_align('center');
    $format{'normal'} ->set_align('vcenter');
    $format{'normal'} ->set_size(12);
    $format{'normal'} ->set_font("Times New Roman");
    $format{'normal'} ->set_border();
    
    $format{'small'} = $workbook->add_format();
    $format{'small'} ->set_align('vcenter');
    $format{'small'} ->set_size(10);
    $format{'small'} ->set_font("Times New Roman");
    $format{'small'} ->set_border();
    
    $format{'seq'} = $workbook->add_format();
    $format{'seq'} ->set_align('vcenter');
    $format{'seq'} ->set_size(11);
    $format{'seq'} ->set_font("Courier New");
    $format{'seq'} ->set_border();
    
    $format{'left'} = $workbook->add_format();
    $format{'left'} ->set_align('vcenter');
    $format{'left'} ->set_size(12);
    $format{'left'} ->set_font("Times New Roman");
    $format{'left'} ->set_border();
    
    $format{'orange'} = $workbook->add_format();
    $format{'orange'} ->set_align('center');
    $format{'orange'} ->set_align('vcenter');
    $format{'orange'} ->set_size(12);
    $format{'orange'} ->set_font("Times New Roman");
    $format{'orange'} ->set_bg_color("#fac090");
    $format{'orange'} ->set_border();

    $format{'skyblue'} = $workbook->add_format();
    $format{'skyblue'} ->set_align('vcenter');
    $format{'skyblue'} ->set_size(12);
    $format{'skyblue'} ->set_font("Times New Roman");
    $format{'skyblue'} ->set_bg_color("#538ed5");
    $format{'skyblue'} ->set_border();

    $format{'bold'} = $workbook->add_format( bold => 1 );
    $format{'blue'} = $workbook->add_format( color => "#538ed5" );
    $format{'redbold'} = $workbook->add_format( color => "#ff0000", bold => 1, );
    $format{'italic'} = $workbook->add_format( italic => 1 );
    $format{'boldblue'} = $workbook->add_format( bold => 1, color => "#538ed5" );
    $format{'bolditalic'} = $workbook->add_format( bold => 1, italic => 1 );
    $format{'blueitalic'} = $workbook->add_format( color => "#538ed5", italic => 1 );
    $format{'boldblueitalic'} = $workbook->add_format( bold => 1, color => "#538ed5", italic => 1 );
    
    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    
    $format{'readme1'} = $workbook->add_format();
    $format{'readme1'}->set_align('center');
    $format{'readme1'}->set_align('vcenter');
    $format{'readme1'}->set_bold();
    $format{'readme1'}->set_size(14);
    $format{'readme1'}->set_font("Times New Roman");
    $format{'readme1'}->set_border();

    $format{'readme2'} = $workbook->add_format();
    $format{'readme2'}->set_align('vcenter');
    $format{'readme2'}->set_bold();
    $format{'readme2'}->set_size(14);
    $format{'readme2'}->set_font("Times New Roman");

    $format{'readme2tmp'} = $workbook->add_format();
    $format{'readme2tmp'}->set_right();

    $format{'readme3'} = $workbook->add_format();
    $format{'readme3'}->set_align('center');
    $format{'readme3'}->set_align('vcenter');
    $format{'readme3'}->set_bold();
    $format{'readme3'}->set_size(11);
    $format{'readme3'}->set_font("Times New Roman");
    $format{'readme3'}->set_border();

    $format{'readme4'} = $workbook->add_format();
    $format{'readme4'}->set_align('vcenter');
    $format{'readme4'}->set_bold();
    $format{'readme4'}->set_size(11);
    $format{'readme4'}->set_font("Times New Roman");
    $format{'readme4'}->set_border();

    $format{'readme5'} = $workbook->add_format();
    $format{'readme5'}->set_align('vcenter');
    $format{'readme5'}->set_size(11);
    $format{'readme5'}->set_font("Times New Roman");
    $format{'readme5'}->set_border();
    $format{'readme5'}->set_text_wrap();

    return %format;
}

# 输出分型结果
sub output_genotype{
    my $hashTyping  = shift @_;
    my $hashMotif   = shift @_;
    my $content     = shift @_;
    my $samples     = shift @_;
    my $output_file = shift @_;

    print LOG "output : $output_file\n";

    my @heads = ('target', @$samples);

    open OUTPUT, ">$output_file";
    print OUTPUT (join "\t", @heads) . "\n";

    foreach my $target_aim_str(sort keys %$hashTyping)
    {

        my @datas;
        foreach my $head(@heads)
        {   
            my $data = "";
            $data = $target_aim_str if($head eq 'target');
            
            if($data eq '' and exists $hashTyping->{$target_aim_str}{$head})
            {
                my @alleles = sort {$a<=>$b} keys %{$hashTyping->{$target_aim_str}{$head}}; # 等位基因
                my @alleles_r; # 等位基因相对量
                my @genotypes; # 分型
                my @genotypes_two_digital; # 2位数表示结果 
                my @alleles_two_digital; # 两位数表示结果
                foreach my $allele(@alleles)
                {   
                    my $allele_two_digital = (length($allele) == 1) ? "0$allele" : $allele; # 客户个性化需求，要求等位基因用两位数字表示，不足的补0
                    map{push @genotypes, $allele; push @genotypes_two_digital, $allele_two_digital} ( 1..$hashTyping->{$target_aim_str}{$head}{$allele}{'Absolute'} );
                    push @alleles_r, $hashTyping->{$target_aim_str}{$head}{$allele}{'Relative'};
                    push @alleles_two_digital, $allele_two_digital;
                }

                # $data = join "/", @alleles;
                # $data = join "/", @alleles_r if($content eq 'allele_r');
                # $data = join "/", @genotypes if($content eq 'genotype');

                # 两位数表示版本
                $data = join "", @alleles_two_digital;
                $data = join "/", @alleles_r if($content eq 'allele_r');
                $data = join "", @genotypes_two_digital if($content eq 'genotype');

            }
            push @datas, $data;
        }
        print OUTPUT (join "\t", @datas) . "\n";
    }
    close OUTPUT;

}

# 输出相对、绝对分型结果
sub output_type_relative_absolute{
    my $hashTyping  = shift @_;
    my $hashMotif   = shift @_;
    my $content     = shift @_;
    my $samples     = shift @_;
    my $output_file = shift @_;

    print LOG "output : $output_file\n";

    # 获取分型结果中，含有的最长motif数量
    my $max_motif_count = 0; 
    foreach my $target_aim_str(keys %$hashTyping)
    {
        foreach my $sample(@$samples)
        {
            foreach my $str_motif_count(keys %{$hashTyping->{$target_aim_str}{$sample}})
            {
                $max_motif_count = $str_motif_count if($str_motif_count > $max_motif_count);
            }
        }
    }

    my @heads = ('Sample', 'target', 1..$max_motif_count);
    open OUTPUT, ">$output_file";
    print OUTPUT (join "\t", @heads) . "\n";

    foreach my $sample(@$samples)
    {
        foreach my $target_aim_str(sort keys %$hashTyping)
        {

            my @datas;
            foreach my $head(@heads)
            {   
                my $data = "";
                $data = $sample                                         if($head eq 'Sample');
                $data = $target_aim_str                                 if($head eq 'target');
                $data = $hashTyping->{$target_aim_str}{$sample}{$head}{$content} if($data eq "" and exists $hashTyping->{$target_aim_str}{$sample} and $hashTyping->{$target_aim_str}{$sample}{$head} );
                push @datas, $data;
            }
            print OUTPUT (join "\t", @datas) . "\n";
        }
    }
    close OUTPUT;
}

# 输出扩增校正比例
sub output_ampr{
    my $hashAMPR      = shift @_;
    my $hashMotif     = shift @_;
    my $output_file   = shift @_;

    print LOG "output : $output_file\n";
 
    my @heads = ('target', 
        'str', 
        'AMPR_mean', 
        'AMPR_count', 
        'AMPR_median', 
        'AMPR_sd', 
        'AMPR_min', 
        'AMPR_max',
        'Detail',
        );
    open OUTPUT, ">$output_file";
    print OUTPUT (join "\t", @heads) . "\n";

    foreach my $target_aim_str(sort keys %$hashAMPR)
    {


        my $motif_seq       = $hashMotif->{$target_aim_str}{'motif'};        # motif序列       
        my $motif_length    = $hashMotif->{$target_aim_str}{'motif_length'}; # motif序列长度

        my @str_identifys     = sort keys %{$hashAMPR->{$target_aim_str}}; # 当前所有STR类型
        my %hashSTRMotifCount = map{ my $motif_count = length($_) / $motif_length; ($_, $motif_count) } @str_identifys; # 每种STR序列含有的motif数量

        foreach my $str_identify(@str_identifys)
        {
            my @datas;
            foreach my $head(@heads)
            {   
                my $data = "";
                $data = $target_aim_str                                 if($head eq 'target');
                $data = "$motif_seq($hashSTRMotifCount{$str_identify})" if($head eq 'str');
                $data = $hashAMPR->{$target_aim_str}{$str_identify}{$head} if($data eq "" and exists $hashAMPR->{$target_aim_str}{$str_identify}{$head} );
                push @datas, $data;
            }
            print OUTPUT (join "\t", @datas) . "\n";
        }
    }
    close OUTPUT;
}

# 输出前滑移比例
sub output_slip_ratio{
    my $hashSlipRatio = shift @_;
    my $hashMotif     = shift @_;
    my $output_file   = shift @_;

    print LOG "output : $output_file\n";

    my @heads = ('target', 
        'str', 
        'SlipRatio', 
        'Type', 
        'SampleCount', 
        'SlipRatio_median', 
        'SlipRatio_sd', 
        'SlipRatio_min', 
        'SlipRatio_max', 
        'Detail',
        );
    open OUTPUT, ">$output_file";
    print OUTPUT (join "\t", @heads) . "\n";

    foreach my $target_aim_str(sort keys %$hashSlipRatio)
    {

        my $motif_seq       = $hashMotif->{$target_aim_str}{'motif'};        # motif序列       
        my $motif_length    = $hashMotif->{$target_aim_str}{'motif_length'}; # motif序列长度

        my @str_identifys     = sort keys %{$hashSlipRatio->{$target_aim_str}}; # 当前所有STR类型
        my %hashSTRMotifCount = map{ my $motif_count = length($_) / $motif_length; ($_, $motif_count) } @str_identifys; # 每种STR序列含有的motif数量

        foreach my $str_identify(@str_identifys)
        {
            my @datas;
            foreach my $head(@heads)
            {   
                my $data = "";
                $data = $target_aim_str                                 if($head eq 'target');
                $data = "$motif_seq($hashSTRMotifCount{$str_identify})" if($head eq 'str');
                $data = $hashSlipRatio->{$target_aim_str}{$str_identify}{$head} if($data eq "" and exists $hashSlipRatio->{$target_aim_str}{$str_identify}{$head} );
                push @datas, $data;
            }
            print OUTPUT (join "\t", @datas) . "\n";
        }
    }
    close OUTPUT;
}


# 输出reads数、频率数据 
sub output_reads_freq{
    my $hashSTR     = shift @_;
    my $hashMotif   = shift @_;
    my $content     = shift @_;
    my $samples     = shift @_;
    my $output_file = shift @_;

    print LOG "output : $output_file\n";

    my @heads = ('target', 'str', @$samples);
    open OUTPUT, ">$output_file";
    print OUTPUT (join "\t", @heads) . "\n";

    foreach my $target_aim_str(sort keys %$hashSTR)
    {

        my $motif_seq       = $hashMotif->{$target_aim_str}{'motif'};        # motif序列       
        my $motif_length    = $hashMotif->{$target_aim_str}{'motif_length'}; # motif序列长度

        my @str_identifys     = sort keys %{$hashSTR->{$target_aim_str}}; # 当前所有STR类型
        my %hashSTRMotifCount = map{ my $motif_count = length($_) / $motif_length; ($_, $motif_count) } @str_identifys; # 每种STR序列含有的motif数量

        foreach my $str_identify(@str_identifys)
        {
            my @datas;
            foreach my $head(@heads)
            {   
                my $data = "";
                $data = $target_aim_str                                 if($head eq 'target');
                $data = "$motif_seq($hashSTRMotifCount{$str_identify})" if($head eq 'str');
                $data = $hashSTR->{$target_aim_str}{$str_identify}{$head}{$content} if($data eq "" and exists $hashSTR->{$target_aim_str}{$str_identify}{$head} and exists $hashSTR->{$target_aim_str}{$str_identify}{$head}{$content} );
                push @datas, $data;
            }
            print OUTPUT (join "\t", @datas) . "\n";
        }
    }
    close OUTPUT;

}

# 分型
sub typing{
    my $hashSTRClean      = shift @_;
    my $hashMotif         = shift @_;
    my $hashAMPR          = shift @_;
    my $samples           = shift @_;
    my $min_relative_copy = shift @_;

    $min_relative_copy = 0 if(not defined $min_relative_copy); # 最低相对拷贝数，默认为0   

    print "5. Typing ... ";  
    my %hashTyping;
    foreach my $target_aim_str(keys %$hashSTRClean)
    {
        my @str_identifys                = sort keys %{$hashSTRClean->{$target_aim_str}}; # 当前所有STR类型

        my $noise_cutoff = $hashMotif->{$target_aim_str}{'noise_cutoff'}; # 设定的噪音阈值
        my $type_cutoff  = $hashMotif->{$target_aim_str}{'type_cutoff'};  # 设定的分型阈值
        my $motif_seq    = $hashMotif->{$target_aim_str}{'motif'};        # motif序列       
        my $motif_length = $hashMotif->{$target_aim_str}{'motif_length'}; # motif序列长度       
        my $ploid        = $hashMotif->{$target_aim_str}{'ploid'};        # 倍体数       
        my $homology     = $hashMotif->{$target_aim_str}{'homology'};     # 同源性

        my $ploidy       = $ploid * $homology; # 当前片段倍体数

        my %hashSTRMotifCount = map{ my $motif_count = length($_) / $motif_length; ($_, $motif_count) } @str_identifys; # 每种STR序列含有的motif数量

        foreach my $sample(@$samples)
        {   
            # (1) 候选分型STR筛查
            my %hashCandidate;
            my $sum = 0;
            foreach my $str_identify(@str_identifys)
            {
                my $freq = exists $hashSTRClean->{$target_aim_str}{$str_identify}{$sample} ? $hashSTRClean->{$target_aim_str}{$str_identify}{$sample}{'Freq'} : 0;
                next if($freq < $noise_cutoff); # 去掉噪音
                $hashCandidate{$str_identify} = $freq;
                $sum += $freq;
            }
            next if($sum == 0); # 没有数据
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum; } keys %hashCandidate; # 重新计算频率

            # （2）复查，去掉滑移
            foreach my $str_candidate(sort keys %hashCandidate){
                next if($hashCandidate{$str_candidate} >= $type_cutoff); # 不用考虑分型线以上的str                
                delete $hashCandidate{$str_candidate} if(exists($hashCandidate{$str_candidate.$motif_seq}) and $hashCandidate{$str_candidate.$motif_seq} > (2 * $type_cutoff)); # 当前STR频率低于分型线，且上一级STR频率高于分型线2倍以上，则当前STR判定为滑移STR，去掉
            }
            $sum = 0;
            map{ $sum += $hashCandidate{$_}; } keys %hashCandidate;
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum } keys %hashCandidate;  

            # （3）去掉多余分型,优先排除频数较低的STR
            my @str_lefts=sort {$hashCandidate{$a}<=>$hashCandidate{$b} or $a cmp $b} keys %hashCandidate;
            if(scalar(@str_lefts) > $ploidy){
               my $remove_count = scalar(@str_lefts) - $ploidy;
               foreach my $remove_index(0..$remove_count-1){
                   delete $hashCandidate{$str_lefts[$remove_index]};
               }
            }
            $sum = 0;
            map{ $sum += $hashCandidate{$_}; } keys %hashCandidate;
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum } keys %hashCandidate; 

            # （4） 扩增百分比校正
            foreach my $str_candidate(sort keys %hashCandidate){
                my $relative_copy  = $hashCandidate{$str_candidate} * $ploidy;
                my $correct_copy   = (exists $hashAMPR->{$target_aim_str}{$str_candidate}) ? $relative_copy / $hashAMPR->{$target_aim_str}{$str_candidate}{'AMPR_mean'} : $relative_copy;
                $hashCandidate{$str_candidate} = $correct_copy; 
            }
            $sum = 0;
            map{ $sum += $hashCandidate{$_}; } keys %hashCandidate;
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum } keys %hashCandidate; 

            # （5）最低相对拷贝数过滤
            foreach my $str_candidate(keys %hashCandidate){
                my $relative_copy = $hashCandidate{$str_candidate} * $ploidy;
                delete $hashCandidate{$str_candidate} if($relative_copy < $min_relative_copy);    
            }
            $sum = 0;
            map{ $sum += $hashCandidate{$_}; } keys %hashCandidate;
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum } keys %hashCandidate;  

            # （5）分型
            my $type_sum = 0;
            my %hashType;
            foreach my $str_candidate(sort keys %hashCandidate){
                my $relative_copy = sprintf "%0.2f", $hashCandidate{$str_candidate} * $ploidy;
                $hashType{$str_candidate}{'Relative'} = $relative_copy; # 相对拷贝数
                $hashType{$str_candidate}{'Absolute'} = ($relative_copy > 1) ? sprintf "%0.0f", $relative_copy : 1; # 绝对拷贝数，<1时设为1，>1时，四舍五入
                $hashType{$str_candidate}{'Differ'}   = ($relative_copy > 1) ? ($hashType{$str_candidate}{'Absolute'} - $hashType{$str_candidate}{'Relative'}) / $hashType{$str_candidate}{'Absolute'} : 0;# 单拷贝造成的差异                
                $type_sum += $hashType{$str_candidate}{'Absolute'};
            }

            # （6）校正
            adjust(\%hashType, $ploidy, $type_sum) if($type_sum != $ploidy); # 分型拷贝数与倍型数不等，校正   

            # （7）数据保存
            foreach my $str_final(keys %hashCandidate){
                my $str_motif_count = $hashSTRMotifCount{$str_final};
                $hashTyping{$target_aim_str}{$sample}{$str_motif_count}{'Relative'} = $hashType{$str_final}{'Relative'};
                $hashTyping{$target_aim_str}{$sample}{$str_motif_count}{'Absolute'} = $hashType{$str_final}{'Absolute'};
            }    
        }
    }

    print "OK\n";
    return %hashTyping;
}

# 校正
sub adjust{
    my $hashType  = shift @_;
    my $ploidy    = shift @_;
    my $type_sum  = shift @_; # 绝对拷贝数和
    my @str_finals = keys %$hashType;

    # 循环，指导拷贝数与倍型数相等
    while($type_sum ne $ploidy){

       if($type_sum < $ploidy){ # 低拷贝,需要增加，优先增加相对拷贝数损失最大的
          my @str_tmps = sort { $hashType->{$a}{'Differ'}   <=> $hashType->{$b}{'Differ'}   } keys %$hashType; # 从小到大排列，理论上，最小值是负值
             @str_tmps = sort { $hashType->{$b}{'Absolute'} <=> $hashType->{$a}{'Absolute'} } keys %$hashType if($hashType->{$str_tmps[0]}{'Differ'} >= 0); # 如果最小值是非负数，则最大拷贝数加1
          my $relative_copy = $hashType->{$str_tmps[0]}{'Relative'};
          my $differ        = $hashType->{$str_tmps[0]}{'Differ'};
          $hashType->{$str_tmps[0]}{'Absolute'} += 1; # 第一个STR增加1拷贝
          $hashType->{$str_tmps[0]}{'Differ'}    = ($relative_copy > 1 ) ? ($hashType->{$str_tmps[0]}{'Absolute'} - $relative_copy) / $hashType->{$str_tmps[0]}{'Absolute'} : ($differ + 1);

       }else{ # 高拷贝
          my @str_tmps = sort { $hashType->{$b}{'Differ'}   <=> $hashType->{$a}{'Differ'}   } keys %$hashType; # 从大到小排列，理论上，最大值是正值
             @str_tmps = sort { $hashType->{$b}{'Absolute'} <=> $hashType->{$a}{'Absolute'} } keys %$hashType if($hashType->{$str_tmps[0]}{'Differ'} <= 0);#如果最大的差异非正数，则最大拷贝数减一
          my $relative_copy = $hashType->{$str_tmps[0]}{'Relative'};
          my $differ        = $hashType->{$str_tmps[0]}{'Differ'};
          $hashType->{$str_tmps[0]}{'Absolute'} -= 1;
          $hashType->{$str_tmps[0]}{'Differ'}    = ($relative_copy > 1) ? ($hashType->{$str_tmps[0]}{'Absolute'} - $relative_copy) / $hashType->{$str_tmps[0]}{'Absolute'} : ($differ - 1);     
       }
       # 重新计算绝对拷贝数和
       $type_sum = 0;
       map{ $type_sum += $hashType->{$_}{'Absolute'}; } @str_finals;
    }   
}

# 计算扩增比例 
sub calculate_amplification_ratio
{
    my $hashSTRClean      = shift @_;
    my $hashMotif         = shift @_;
    my $samples           = shift @_;
    my $min_relative_copy = shift @_;

    $min_relative_copy = 0 if(not defined $min_relative_copy); # 最低相对拷贝数，默认为0 

    print "4. Calculate amplification ratio ... ";
    
    # 一：简单分型，学习
    my %hashLearn;
    foreach my $target_aim_str(keys %$hashSTRClean)
    {
        my @samples_depth_ok             = get_depth_enough_sample($hashSTRClean->{$target_aim_str}, $samples); # 获取深度合适的样本
        my @str_identifys                = sort keys %{$hashSTRClean->{$target_aim_str}}; # 当前所有STR类型

        my $noise_cutoff = $hashMotif->{$target_aim_str}{'noise_cutoff'}; # 设定的噪音阈值
        my $type_cutoff  = $hashMotif->{$target_aim_str}{'type_cutoff'};  # 设定的分型阈值
        my $motif_seq    = $hashMotif->{$target_aim_str}{'motif'};        # motif序列       
        my $ploid        = $hashMotif->{$target_aim_str}{'ploid'};        # 倍体数       
        my $homology     = $hashMotif->{$target_aim_str}{'homology'};     # 同源性

        my $ploidy       = $ploid * $homology; # 当前片段倍体数

        foreach my $sample(@samples_depth_ok)
        {   
            # (1) 候选分型STR筛查
            my %hashCandidate;
            my $sum = 0;
            foreach my $str_identify(@str_identifys)
            {
                my $freq = exists $hashSTRClean->{$target_aim_str}{$str_identify}{$sample} ? $hashSTRClean->{$target_aim_str}{$str_identify}{$sample}{'Freq'} : 0;
                next if($freq < $noise_cutoff); # 去掉噪音
                $hashCandidate{$str_identify} = $freq;
                $sum += $freq;
            }
            next if($sum == 0); # 没有数据
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum; } keys %hashCandidate; # 重新计算频率

            # （2）复查，去掉滑移
            foreach my $str_candidate(sort keys %hashCandidate){
                next if($hashCandidate{$str_candidate} >= $type_cutoff); # 不用考虑分型线以上的str                
                delete $hashCandidate{$str_candidate} if(exists($hashCandidate{$str_candidate.$motif_seq}) and $hashCandidate{$str_candidate.$motif_seq} > (2 * $type_cutoff)); # 当前STR频率低于分型线，且上一级STR频率高于分型线2倍以上，则当前STR判定为滑移STR，去掉
            }
            $sum = 0;
            map{ $sum += $hashCandidate{$_}; } keys %hashCandidate;
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum } keys %hashCandidate;  

            # （3）检查剩余的STR数量
            my @str_lefts = keys %hashCandidate;
            next if(scalar(@str_lefts) > $ploidy);  # 超过倍体数量，跳过

            # （4）最低相对拷贝数过滤
            foreach my $str_candidate(keys %hashCandidate){
                my $relative_copy = $hashCandidate{$str_candidate} * $ploidy;
                delete $hashCandidate{$str_candidate} if($relative_copy < $min_relative_copy);    
            }
            $sum = 0;
            map{ $sum += $hashCandidate{$_}; } keys %hashCandidate;
            map{ $hashCandidate{$_} = $hashCandidate{$_} / $sum } keys %hashCandidate;  

            # （5）分型
            my $type_sum = 0;
            my %hashType;
            foreach my $str_candidate(sort keys %hashCandidate){
                my $relative_copy = sprintf "%0.2f", $hashCandidate{$str_candidate} * $ploidy;
                $hashType{$str_candidate}{'Relative'} = $relative_copy; # 相对拷贝数
                $hashType{$str_candidate}{'Absolute'} = ($relative_copy > 1) ? sprintf "%0.0f", $relative_copy : 1; # 绝对拷贝数，<1时设为1，>1时，四舍五入
                $type_sum += $hashType{$str_candidate}{'Absolute'};
            }
            next if($type_sum != $ploidy);# 过滤掉不能学习校正的样本   

            # （6）记录学习用数据
            foreach my $str_final(keys %hashType){
                my $copy1_perc = sprintf "%0.5f", $hashType{$str_final}{'Relative'} / $hashType{$str_final}{'Absolute'}; # 一拷贝所占相对拷贝数
                $hashLearn{$target_aim_str}{$str_final}{$sample} = $copy1_perc;
            }   
        }

    }

    # 二：校正计算  
    my %hashAMPR; # 扩增校正，单拷贝比例, amplication ratio
    foreach my $target_aim_str(keys %hashLearn){
        foreach my $str_final(sort keys %{$hashLearn{$target_aim_str}}){

            # (1) 计算离群阈值
            my @percs = map{ $hashLearn{$target_aim_str}{$str_final}{$_}; } keys %{$hashLearn{$target_aim_str}{$str_final}};

            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@percs);
            my $quan25 = (defined $stat->quantile(1) ) ? $stat->quantile(1) : "NA"; # 下4分位
            my $quan75 = (defined $stat->quantile(3) ) ? $stat->quantile(3) : "NA"; # 上4分位
            my $down_cutoff = $quan25 - 1.5 * ($quan75 - $quan25); # 离群上下阈值
            my $up_cutoff   = $quan75 + 1.5 * ($quan75 - $quan25);

            # （2）去掉离群数据
            my @no_outline_infos;
            my @no_outline_percs;
            foreach my $sample(keys %{$hashLearn{$target_aim_str}{$str_final}})
            {
                next if($hashLearn{$target_aim_str}{$str_final}{$sample} < $down_cutoff or $hashLearn{$target_aim_str}{$str_final}{$sample} > $up_cutoff);                
                my $perc = $hashLearn{$target_aim_str}{$str_final}{$sample};
                push @no_outline_infos, "$sample=$perc";
                push @no_outline_percs, $perc;
            }

            # （3）统计
            my $stat2 = Statistics::Descriptive::Full->new();
            $stat2->add_data(@no_outline_percs);
            $hashAMPR{$target_aim_str}{$str_final}{'AMPR_mean'}   = $stat2->mean();
            $hashAMPR{$target_aim_str}{$str_final}{'AMPR_count'}  = $stat2->count();
            $hashAMPR{$target_aim_str}{$str_final}{'AMPR_median'} = $stat2->median();
            $hashAMPR{$target_aim_str}{$str_final}{'AMPR_sd'}     = $stat2->standard_deviation();
            $hashAMPR{$target_aim_str}{$str_final}{'AMPR_min'}    = $stat2->min();
            $hashAMPR{$target_aim_str}{$str_final}{'AMPR_max'}    = $stat2->max();
            $hashAMPR{$target_aim_str}{$str_final}{'Detail'}      = join ";", @no_outline_infos; # 样本详情
        }
    }

    print "OK\n";
    return %hashAMPR;
}

# 滑移校正
sub slip_correction{
    my $hashSTR       = shift @_;
    my $hashMotif     = shift @_;
    my $hashSlipRatio = shift @_;
    my $samples       = shift @_;

    print "3. Slip correction ... ";
    # (1) 计算每个样本的每种STR应该去掉多少reads
    my %hashSTRMinus;
    foreach my $target_aim_str(keys %$hashSTR)
    {
        my $motif_seq                    = $hashMotif->{$target_aim_str}{'motif'};        # motif序列

        foreach my $str_identify(keys %{$hashSTR->{$target_aim_str}})
        {   
            my $str_identify_add_1motif = "$str_identify$motif_seq"; # 多一个motif的STR序列
            my $slip_ratio              = (exists $hashSlipRatio->{$target_aim_str}{$str_identify_add_1motif}) ? $hashSlipRatio->{$target_aim_str}{$str_identify_add_1motif}{'SlipRatio'} : 0; # 滑移比例
            foreach my $sample(@$samples)
            {   
                # 样本在当前与多一个motif的STR序列上都有数据
                next if(not exists $hashSTR->{$target_aim_str}{$str_identify}{$sample});
                next if(not exists $hashSTR->{$target_aim_str}{$str_identify_add_1motif} or not exists $hashSTR->{$target_aim_str}{$str_identify_add_1motif}{$sample});

                $hashSTRMinus{$target_aim_str}{$str_identify}{$sample} = sprintf "%0.0f", $slip_ratio * $hashSTR->{$target_aim_str}{$str_identify_add_1motif}{$sample}{'Reads'}; # 当前str获得的上一级滑移的reads数量 
            }
        }
    }

    # (2) 去掉滑移reads数
    my %hashSTRClean;
    my %hashTargetSum;
    foreach my $target_aim_str(keys %$hashSTR)
    {
        foreach my $str_identify(keys %{$hashSTR->{$target_aim_str}})
        {   
            foreach my $sample(@$samples)
            {   
                # 样本在当前与多一个motif的STR序列上都有数据
                next if(not exists $hashSTR->{$target_aim_str}{$str_identify}{$sample});
                my $minus = exists $hashSTRMinus{$target_aim_str}{$str_identify}{$sample} ? $hashSTRMinus{$target_aim_str}{$str_identify}{$sample} : 0; # 最长的STR是无法进行前滑移校正的
                my $reads_count_clean = $hashSTR->{$target_aim_str}{$str_identify}{$sample}{'Reads'} - $minus;
                $reads_count_clean    = 0 if($reads_count_clean < 0); # 可能会是负数

                $hashSTRClean{$target_aim_str}{$str_identify}{$sample}{'Reads'}  = $reads_count_clean; # 当前str获得的上一级滑移的reads数量 
                $hashTargetSum{$target_aim_str}{$sample}                        += $reads_count_clean; # 样本在片段上的reads总数
            }
        }
    }  
    # (3) 重新计算频率
    foreach my $target_aim_str(keys %$hashSTR)
    {
        foreach my $str_identify(keys %{$hashSTR->{$target_aim_str}})
        {   
            foreach my $sample(@$samples)
            {   
                # 样本在当前与多一个motif的STR序列上都有数据
                next if(not exists $hashSTR->{$target_aim_str}{$str_identify}{$sample});
                my $freq = sprintf "%0.5f", $hashSTRClean{$target_aim_str}{$str_identify}{$sample}{'Reads'} / $hashTargetSum{$target_aim_str}{$sample};
                   $freq = 0 if($hashSTR->{$target_aim_str}{$str_identify}{$sample}{'Freq'} == 0); # 深度不足

                $hashSTRClean{$target_aim_str}{$str_identify}{$sample}{'Freq'} = $freq; 
            }
        }
    }    

    print "OK\n";
    return %hashSTRClean;
}

# 滑移比例计算
sub calculate_slip_ratio{
    my $hashSTR   = shift @_;
    my $hashMotif = shift @_;
    my $samples   = shift @_;

    print "2. Calculate Slip Ratio ... ";
    # （1）获得滑移比例候选数据组
    my %hashSlipCandidate; 
    foreach my $target_aim_str(sort keys %$hashSTR){
        my @samples_depth_ok             = get_depth_enough_sample($hashSTR->{$target_aim_str}, $samples); # 获取深度合适的样本
        my @str_identifys                = sort keys %{$hashSTR->{$target_aim_str}}; # 当前所有STR类型

        my $noise_cutoff = $hashMotif->{$target_aim_str}{'noise_cutoff'}; # 设定的噪音阈值
        my $type_cutoff  = $hashMotif->{$target_aim_str}{'type_cutoff'};  # 设定的分型阈值
        my $motif_seq    = $hashMotif->{$target_aim_str}{'motif'};        # motif序列
        
        foreach my $sample(@samples_depth_ok)
        {
            foreach my $str_identify(@str_identifys)
            {   
                # 当前样本当前str序列频率满足条件
                next if(not exists $hashSTR->{$target_aim_str}{$str_identify}{$sample});
                my $freq = $hashSTR->{$target_aim_str}{$str_identify}{$sample}{'Freq'};
                next if($freq < $type_cutoff); # 频率低于分型阈值，不考虑

                # 当前样本当前str前移一个motif的序列频率满足条件
                my $front_str_identify = get_front_str($str_identify, $motif_seq); # 获取当前STR序列少一个motif的序列
                next if($front_str_identify !~ /\w/); # 跳过1motif 
                my $front_str_identify_freq = (exists $hashSTR->{$target_aim_str}{$front_str_identify} and $hashSTR->{$target_aim_str}{$front_str_identify}{$sample} ) ? $hashSTR->{$target_aim_str}{$front_str_identify}{$sample}{'Freq'} : 0; # 前移str的频率
                next if($front_str_identify_freq == 0 or $front_str_identify_freq > $noise_cutoff); # 前移STR高于噪音阈值

                # 滑移比例计算
                $hashSlipCandidate{$target_aim_str}{$str_identify}{$sample} = sprintf "%0.5f", $front_str_identify_freq / $freq;
            }
        }      
    }

    # （2）计算滑移比例均值，去掉离群值
    my %hashSlipRatio;
    foreach my $target_aim_str(keys %hashSlipCandidate)
    {
        foreach my $str_identify(keys %{$hashSlipCandidate{$target_aim_str}})
        {
            # (1) 计算离群阈值
            my @ratios = map{ $hashSlipCandidate{$target_aim_str}{$str_identify}{$_} } keys %{$hashSlipCandidate{$target_aim_str}{$str_identify}};
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@ratios);
            my $quan25 = (defined $stat->quantile(1) ) ? $stat->quantile(1) : "NA"; # 下4分位
            my $quan75 = (defined $stat->quantile(3) ) ? $stat->quantile(3) : "NA"; # 上4分位
            my $down_cutoff = $quan25 - 1.5 * ($quan75 - $quan25); # 离群上下阈值
            my $up_cutoff   = $quan75 + 1.5 * ($quan75 - $quan25);

            # （2）去掉离群数据
            my @no_outline_infos;
            my @no_outline_percs;
            foreach my $sample(keys %{$hashSlipCandidate{$target_aim_str}{$str_identify}})
            {
                next if($hashSlipCandidate{$target_aim_str}{$str_identify}{$sample} < $down_cutoff or $hashSlipCandidate{$target_aim_str}{$str_identify}{$sample} > $up_cutoff);                
                my $ratio = $hashSlipCandidate{$target_aim_str}{$str_identify}{$sample};
                push @no_outline_infos, "$sample=$ratio";
                push @no_outline_percs, $ratio;
            }

            # （3）统计
            my $stat2 = Statistics::Descriptive::Full->new();
            $stat2->add_data(@no_outline_percs);
            $hashSlipRatio{$target_aim_str}{$str_identify}{'Type'}             = 'SampleGenerate';# 数据来源于样本实测
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio'}        = $stat2->mean(); # 滑移比例
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SampleCount'}      = $stat2->count(); # 样本数量
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio_median'} = $stat2->median(); # 中值
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio_sd'}     = $stat2->standard_deviation(); # 标准差
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio_min'}    = $stat2->min(); # 最小值
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio_max'}    = $stat2->max(); # 最大值
            $hashSlipRatio{$target_aim_str}{$str_identify}{'Detail'}           = join ";", @no_outline_infos; # 样本详情         
        }
    }

    # (3) 一元二次方程，拟合滑移比例，补充无法获得滑移比例的STR
    foreach my $target_aim_str(sort keys %$hashSTR){
        my @str_identifys                = sort keys %{$hashSTR->{$target_aim_str}}; # 当前所有STR类型

        my $motif_seq    = $hashMotif->{$target_aim_str}{'motif'};        # motif序列
        my $motif_length = $hashMotif->{$target_aim_str}{'motif_length'}; # motif长度

        # 获取用于拟合的数据
        my @motif_counts = (0); # motif数量
        my @ratios       = (0); # 滑移比例
        foreach my $str_identify(@str_identifys)
        {
            next if(not exists $hashSlipRatio{$target_aim_str}{$str_identify});
            my $str_identify_motif_count = length($str_identify) / $motif_length;
            push @motif_counts, $str_identify_motif_count;
            push @ratios, $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio'};
        }
        # 一元二次拟合，公式为 y = a*x^2 + b*x + c 
        my ($a, $b, $c) = unitary_quadratic_equation(\@motif_counts, \@ratios);
        # 补充缺失值
        foreach my $str_identify(@str_identifys){
            next if(exists $hashSlipRatio{$target_aim_str}{$str_identify});
            my $str_identify_motif_count = length($str_identify) / $motif_length;           

            my $y = $a * $str_identify_motif_count * $str_identify_motif_count + $b * $str_identify_motif_count + $c;# 拟合值
            $hashSlipRatio{$target_aim_str}{$str_identify}{'Type'}        = "unitary_quadratic_equation";#强制为0，有可能会出现负数
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SlipRatio'}   = ($y >= 0) ? $y : 0;#强制为0，有可能会出现负数
            $hashSlipRatio{$target_aim_str}{$str_identify}{'SampleCount'} = ""; # 样本数量
            $hashSlipRatio{$target_aim_str}{$str_identify}{'Detail'}      = ""; # 样本详情 
        }
    }
    print "OK\n";
    return %hashSlipRatio;
}

# 一元二次方程拟合
sub unitary_quadratic_equation{
    my $motif_counts = shift @_;
    my $ratios       = shift @_;
    
    my ($a, $b, $c) = (0, 0, 0);
    return ($a, $b, $c) if(@$motif_counts < 2); # 拟合最低要有两个点
    my $x = join ",",@$motif_counts;
    my $y = join ",",@$ratios;
    my $R = Statistics::R->new();
    $R->startR;
    $R->send(qq` x=c($x)\n`);   
    $R->send(qq` y=c($y)\n`);   
    $R->send(qq` data = data.frame(x,y)\n`);
    $R->send(qq` result = lm(y ~ I(x^2)+x , data = data)\n`);
    $a = $R->get('result$coefficients[2]');# 获取三个系数
    $b = $R->get('result$coefficients[3]');
    $c = $R->get('result$coefficients[1]');
    $R->stop();
    $a = 0 if($a !~ /\d/);
    $b = 0 if($b !~ /\d/);
    $c = 0 if($c !~ /\d/);
    return ($a, $b, $c);
}

# 获取当前STR序列少一个motif的序列
sub get_front_str{
    my $str_identify = shift @_;
    my $motif_seq    = shift @_;

    my $motif_seq_length         = length($motif_seq);
    my $str_identify_length      = length($str_identify);
    my $str_identify_motif_count = $str_identify_length / $motif_seq_length; # str序列含有的motif数量

    my $front_str_identify = ($str_identify_motif_count == 1) ? "" : $motif_seq x ($str_identify_motif_count - 1);# 前移STR;

    return $front_str_identify;
}

# 获取测序深度足够高的样本
sub get_depth_enough_sample
{
    my $hashTargetAimSTR = shift @_;
    my $samples          = shift @_;

    # 当前所有STR类型
    my @str_identifys = keys %$hashTargetAimSTR; 
    # 每个样本的深度
    my %hashSampleDepth; 
    my $sample_count = 0;
    my $sample_depth_sum = 0;
    foreach my $sample(@$samples)
    {   
        my $freq_sum = 0;
        foreach my $str_identify(@str_identifys)
        {
            next if(not exists $hashTargetAimSTR->{$str_identify}{$sample});
            $hashSampleDepth{$sample} += $hashTargetAimSTR->{$str_identify}{$sample}{'Reads'};
            $freq_sum += $hashTargetAimSTR->{$str_identify}{$sample}{'Freq'};
        }
        if($freq_sum == 0) # 去掉不考虑的样本
        {
            delete $hashSampleDepth{$sample};
            next;
        }

        $sample_count++;
        $sample_depth_sum += $hashSampleDepth{$sample};
    }

    # 样本深度均值的50%
    my $depth_mean_50 = ($sample_count > 0 ) ? ($sample_depth_sum / $sample_count) * 0.5 : 0;
    # 合格样本
    my @samples_depth_ok = ();
    map{ push @samples_depth_ok, $_ if($hashSampleDepth{$_} >= $depth_mean_50); } (keys %hashSampleDepth);

    return @samples_depth_ok;
}

# 读入motif信息文件
sub read_motif{
    my $motif_file = shift @_;
    my $hashMotif  = shift @_;

    print "[Reading] $motif_file ... ";
    open MOTIF, $motif_file or die "Error\n[Error] file lost：$motif_file\n";
    my $line1 = <MOTIF>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;

    # 关键表头检测
    my @key_heads = ('target', 'motif', 'ploid', 'homology', 'noise_cutoff', 'type_cutoff');  
    my @lost_heads = grep{  not $_ ~~ @heads} @key_heads;
    die "Error\n[Error] motif file lost key head: @lost_heads\n" if(@lost_heads > 0);
    
    my $row = 1;
    while(<MOTIF>)
    {
        $_=~s/[\r\n]//g;
        $row++;
        next if($_!~/\w/);
        my @datas = split /\t/, $_;

        my %hashTmp;
        foreach my $col(0..$#heads)
        {   
            my $head = $heads[$col];
            my $value = exists $datas[$col] ? $datas[$col] : "";
            $value=~s/^\s+//;
            $value=~s/\s+$//;
            die "Error\n[Error] target name should be empty: $value\nerror line：$row\n"                      if($head eq 'target'       and $value !~ /\w/);
            die "Error\n[Error] Motif must be composed of ATCG nucleotide, and be upcase ：$value\nerror line：$row\n" if($head eq 'motif'        and ($value !~ /\w/ or $value =~ /^ATCG/));
            die "Error\n[Error] ploid must be integer ：$value\nerror line：$row\n"                   if($head eq 'ploid'        and ($value !~ /\w/ or $value !~ /^\d+$/));
            die "Error\n[Error] homology must be integer：$value\nerror line：$row\n"                if($head eq 'homology'     and ($value !~ /\w/ or $value !~ /^\d+$/));
            die "Error\n[Error] noise_cutoff can only be within (0,1) ：$value\nerror line：$row\n" if($head eq 'noise_cutoff'  and ($value !~ /\w/ or $value !~ /^\d+\.\d+$/));
            die "Error\n[Error] type_cutoff can only be within ：$value\nerror line：$row\n"  if($head eq 'type_cutoff'   and ($value !~ /\w/ or $value !~ /^\d+\.\d+$/));
            $hashTmp{$head} = $value;
        }

        my $target = $hashTmp{'target'};
        die "Error\n[Error] noise cutoff must be less than type cutoff ：noise_cutoff = $hashTmp{'noise_cutoff'}, type_cutoff = $hashTmp{'type_cutoff'}\nerror line：$row\n" if($hashTmp{'noise_cutoff'} >= $hashTmp{'type_cutoff'});
        die "Error\n[Error] duplicate target name, which should be unique：$target\nerror line：$row\n" if(exists $hashMotif->{$target});
        print "[Warnings] copy number is over 6, which could lead to error result ：$target\nwarning line：$row\n" if(($hashTmp{'ploid'} * $hashTmp{'homology'}) > 6);

        # 结果保留
        map{ $hashMotif->{$target}{$_} = $hashTmp{$_} } keys %hashTmp;
        $hashMotif->{$target}{'motif_length'} = length($hashTmp{'motif'});
        
    }
    close MOTIF;
    print "OK\n";
}

# col格式文件读入
sub read_str_col{
    my $input     = shift @_;
    my $hashSTR   = shift @_;
    my $hashMotif = shift @_;
    my $samples   = shift @_;
    my $min_reads = shift @_;
    
    print "[Reading] $input [format: col] ... ";
    open INPUT, $input or die "Error\n[Error] lost file：$input\n";
    my $line1 = <INPUT>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    my %hashSTRDepth;  # 记录每个片段上，每个样本的总reads数量
    
    my $row = 1;
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        $row++;
        next if($_!~/\w/);
        my ($target, $str, $sample, $reads_count, $tmp) = split /\t/, $_;
        $target =~ s/^\s+//g;  # 去掉开头、结尾的空格,防止客户失误
        $target =~ s/\s+$//g;  #  
        die "Error\n[Error] target name lost\nerror line：$row\n" if($target !~ /\w/);
 
        die "Error\n[Error] format 'col' contain 4 column at least：target name、STR type、sample name、reads count\n error line：$row\n" if(defined $tmp);

        # 客户定义的str格式不对
        die "Error\n[Error] The second column format  is error in the input file, which must meet format “[ATCG]+\\d+” ，example “TCCA(27)”。\nerror detail：$target\t$str\nerror line：$row\n" if($str !~ /[ATCG]+\(\d+\)/); 

        my ($motif, $motif_count) = $str =~ /([ATCG]+)\((\d+)\)/i;  # 获取motif序列、重复次数
        $motif = uc($motif);
        my $str_seq = $motif x $motif_count;  # 转成str序列

        # 片段没有声明motif信息，或motif序列不一致
        die "Error\n[Error] target lost motif information, please update motif file and reupload \nerror detail：$target\nerror line：$row\n" if(not exists $hashMotif->{$target});
        die "Error\n[Error] target motif sequence is different frome motif file, please recheck and reupload。\nerror detail：$target\t$str\nerror line：$row\n" if($hashMotif->{$target}{'motif'} ne $motif);

        $reads_count =~ s/^\s+//g;  # 去掉开头、结尾的空格
        $reads_count =~ s/\s+$//g;  
        # 非数字
        die "Error\n[Error] non-numeric info is included in file。\nerror detail：$target\t$str\t$sample\t$reads_count\nerror line：$row\n" if($reads_count !~ /^\d+$/);

        $hashSTR->{$target}{$str_seq}{$sample}{'Reads'} = $reads_count;
        $hashSTRDepth{$target}{$sample} += $reads_count;
    }

    #  统计每个样本每种STR的频率
    foreach my $target(keys %$hashSTR)
    {
        foreach my $str_seq(keys %{$hashSTR->{$target}})
        {
            foreach my $sample(@samples)
            {
                next if(not exists $hashSTR->{$target}{$str_seq}{$sample});  # 样本没有这种str序列

                my $perc = sprintf "%0.5f", $hashSTR->{$target}{$str_seq}{$sample}{'Reads'} / $hashSTRDepth{$target}{$sample};
                $perc = 0 if($hashSTRDepth{$target}{$sample} < $min_reads); # 深度不足，频率设为0

                $hashSTR->{$target}{$str_seq}{$sample}{"Freq"} = $perc;
            }
        }
    }
    print "OK\n";
}


# row格式文件读入
sub read_str_row{
    my $input     = shift @_;
    my $hashSTR   = shift @_;
    my $hashMotif = shift @_;
    my $samples   = shift @_;
    my $min_reads = shift @_;
    
    print "[Reading] $input [format: row] ... ";
    open INPUT, $input or die "Error\n[Error] file lost：$input\n";
    my $line1 = <INPUT>;
       $line1 =~ s/[\r\n]//g;
    my @heads = split /\t/, $line1;
    @$samples = (@heads[2..$#heads]);  # 样本
    my %hashSTRDepth;  # 记录每个片段上，每个样本的总reads数量

    my $row = 1;
    while(<INPUT>)
    {
        $_=~s/[\r\n]//g;
        $row++;
        next if($_!~/\w/);
        my @datas = split /\t/, $_;
        my ($target, $str) = ($datas[0], $datas[1]);
        $target =~ s/^\s+//g;  # 去掉开头、结尾的空格,防止客户失误
        $target =~ s/\s+$//g;  #        
        die "Error\n[Error] target name is lost\n error line：$row\n" if($target !~ /\w/);
        die "Error\n[Error] input file error, which should contain 3 columns at least\nerror line：$row" if(not exists $datas[2]);

        # 客户定义的str格式不对
        die "Error\n[Error] The second column format  is error in the input file, which must meet format: “[ATCG]+\\d+” ，example: “TCCA(27)”。\nerror detail：$target\t$str\nerror line：$row\n" if($str !~ /[ATCG]+\(\d+\)/);

        my ($motif, $motif_count) = $str =~ /([ATCG]+)\((\d+)\)/i;  # 获取motif序列、重复次数
        $motif = uc($motif);
        my $str_seq = $motif x $motif_count;  # 转成str序列

        # 片段没有声明motif信息，或motif序列不一致
        die "Error\n[Error] target lost motif infomation, please update motif file and reupload\nerror detail：$target\nerror line：$row\n" if(not exists $hashMotif->{$target});
        die "Error\n[Error] target motif sequence is different frome motif file, please recheck and reupload。\nerror detail：$target\t$str\nmotif信息文件中的motif序列为：$hashMotif->{$target}{'motif'}\nerror line：$row\n" if($hashMotif->{$target}{'motif'} ne $motif);

        # 读入每个样本的reads数量
        foreach my $col(2..$#heads)
        {
            my $sample = $heads[$col];
            my $reads_count = (exists $datas[$col]) ? $datas[$col] : "";
               $reads_count =~ s/^\s+//g;  # 去掉开头、结尾的空格
               $reads_count =~ s/\s+$//g;  #
            next if($reads_count eq "");  # 空白
            # 非数字
            die "Error\n[Error] non-numeric info is included in file。\nerror detail: $target\t$str\t$sample\t$reads_count\nerror line：$row\n" if($reads_count !~ /^\d+$/);

            $hashSTR->{$target}{$str_seq}{$sample}{'Reads'} = $reads_count;
            $hashSTRDepth{$target}{$sample} += $reads_count;
        }
    }

    #  统计每个样本每种STR的频率
    foreach my $target(keys %$hashSTR)
    {
        foreach my $str_seq(keys %{$hashSTR->{$target}})
        {
            foreach my $sample(@samples)
            {
                next if(not exists $hashSTR->{$target}{$str_seq}{$sample});  # 样本没有这种str序列

                my $perc = sprintf "%0.5f", $hashSTR->{$target}{$str_seq}{$sample}{'Reads'} / $hashSTRDepth{$target}{$sample};
                $perc = 0 if($hashSTRDepth{$target}{$sample} < $min_reads); # 深度不足，频率设为0

                $hashSTR->{$target}{$str_seq}{$sample}{"Freq"} = $perc;
            }
        }
    }
    print "OK\n";
}

sub help{
    my $info = "
Program: str type， A typing tool based on NGS SSR abundance table
Version: 2019-12-17
Contact: 129 bin gan

Usage:   perl ".(File::Spec->splitpath(File::Spec->rel2abs($0)))[2]." [options]

Options:
        [required]
         --input/-i         Sample STR abundance table（Supports simultaneous analysis of multiple fragments）
                            you can get it in the 'str_count.txt' file which produced by my software 'SSRseq_count' (you also need to do a little bit change to fit str_type.pl software)
         --motif/-m         STR motif information table
                            this file declares each fragment's motif information, copy number information, etc. (must include all fragments appearing in input, the header is fixed)
                            Must contain 6 columns, the column names are:
                            target：fragment name (please keep it same as the Sample STR abundance table )
                            motif：The minimum constituent unit of STR fragment. It is a short sequence composed of capitalized ATCG which must be consistent with the motif in the input file.
                            ploid：ploidy number of species
                            homology：The homology number of this fragment in the genome, usually set to 1. If there are n homologs (there are n different positions in the genome have exactly the same sequence), then set n
                                    Note: copy number of final type = ploid * homology
                            noise_cutoff：During copy number analysis, the noise threshold (STRs with frequencies below the this threshold will be directly excluded), usually set to 0.6 * type_cutoff
                            type_cutoff：During copy number analysis, the typing threshold (STRs with frequencies higher than this threshold will be considered true), usually set to 0.5 * (1 / (ploid * homology))
                                        The STR between the noise threshold and the typing threshold will be corrected and typed by a series of algorithms.

         --format/-f        the format of Sample STR abundance table：row or col
                            row format：each row records the number of reads of STR on a certain fragment for every sample. If a sample contains no reads, it can be 0 or left blank
                            First column: fragment name
                            Second column: STR sequence type
                            The third column and after: the number of reads of this STR sequence contained in each sample
                            The file must include a table header. You can give the first two columns any names, but we recommend 'target' and 'str'. The following columns should be named by sample names.
                            example file： input_example_row.txt

                            col format：each line records the number of reads of STR on a certain fragment for a certain sample. If the sample contains no reads, it can be 0 or this line can be deleted, blank is not allowed. Each line must have 4 columns.
                            First column: fragment name
                            Second column: STR sequence type
                            Third column: sample name
                            Fourth column: number of reads
                            example file： input_example_col.txt

                            Matters needing attention：
                            （a）Fragment name: choose any name, but use letters, numbers, or underscores only, spaces is not allowed. The fragment name must be unique (the STR of the same target fragment must have the same fragment name), otherwise unexpected errors will occur.
                            （b）STR sequence type: Genesky uses a fixed format expressed as 'motif (n)', which means this STR is a sequence composed of n motif. motif is a short sequence of ATCG and must be capitalized.e.g. AGT(8)
         --output_dir/-o     output directory

        [optional]
        --min_reads          the minimum reads required for typing, default: 30
        --min_rc             the minimum relative copy number required for typing, ranges from 0 to 1. if the relative copy number of the sample is lower than this value, alleles are directly excluded, default: 0

         --help/-h           help doc
    \n";
    return $info;
}

