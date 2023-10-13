---
layout:     post
title:      使用googlecloudsdk下载GATK数据
date:       2023-10-13
author:     champeil
catalog: true
tags:
    - software
    - googlecloud
    - download
    - gatk
    - reference
    - window
---
# 前言
- 在[gatk bundle](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=true)中，gatk使用googlecloud进行数据储存，在网页界面上只能单一文件进行下载
- Google Cloud CLI可以在本地对储存在gs中的数据进行批量下载

# 使用
## 安装Google Cloud CLI
- 下载客户端（windows）[Google Cloud CLI](https://dl.google.com/dl/cloudsdk/channels/rapid/GoogleCloudSDKInstaller.exe?hl=zh-cn)
- 安装
- 打开Google Cloud SDK Shell
- 初始化并配置代理（clash）
- 使用gstuil下载（在bundle中多选文件，点击下载，就会有一个gstuil命令显示，直接将命令复制到shell中运行）

```shell
gcloud init # set the init information [configure info, google account varification]
# error would occur cause the firewall, so need to set proxy
# error: ERROR: gcloud crashed (ConnectionError): HTTPSConnectionPool(host='oauth2.googleapis.com', port=443): Max retries exceeded with url: /token (Caused by NewConnectionError('<urllib3.connection.HTTPSConnection object at 0x00000198B0B26970>: Failed to establish a new connection: [WinError 10060] A connection attempt failed because the connected party did not properly respond after a period of time, or established connection failed because connected host has failed to respond'))

# reference: https://cloud.google.com/sdk/docs/proxy-settings?hl=zh-cn
gcloud config set proxy/type [PROXY_TYPE][http,https,socks4,socks5]
gcloud config set proxy/address [PROXY_IP_ADDRESS][your local ip if you connected the proxy]
gcloud config set proxy/port [PROXY_PORT][the port of the proxy]

# then init again
gcloud init

# download example 
gsutil -m cp "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz" .

```

